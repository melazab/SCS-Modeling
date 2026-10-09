"""Tetrahedron face adjacency, computed once per mesh in memory.

face_neighbours() returns, for every face of every tet, the tet on the other
side (-1 on the mesh boundary). Everything that used to re-sort faces -- the
dura check, the preview's tissue interfaces, Elmer's boundary export -- reads
this one array instead.

WHY NOT THE OLD DISK-PARTITIONED SORT (replaced 2026-10-09)
The previous face_batches() hashed faces into ~4 x tets / 1e6 bucket files and
appended to them 100 000 tets at a time: ~435 000 small appends for a 52 M-tet
mesh. On Pioneer's NFS home one append costs ~16 ms (0.07 ms on local disk),
so each sort took ~1 h and the fine dorsal run did four of them. Here faces
are bucketed in memory by their smallest node id and each bucket is sorted on
one packed 64-bit key, in threads (numpy releases the GIL while sorting).
Memory is ~110 bytes per tet: ~6 GB at 52 M tets, ~2 GB at 17 M.
"""
from concurrent.futures import ThreadPoolExecutor

import numpy as np

FACE_NODES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
_FACE_NODES = np.array(FACE_NODES)
# Smallest-node-id range per bucket: ~1.5 M faces per bucket at typical density.
BUCKET_BITS = 16


def face_neighbours(tets, threads=1):
    """(T, 4) array: the tet across face FACE_NODES[k] of tet t, or -1."""
    tets = np.asarray(tets)
    count = len(tets)
    nbr = np.full((count, 4), -1, dtype=np.int32 if 4 * count < 2**31 else np.int64)
    if count == 0:
        return nbr
    bits = max(1, int(tets.max()).bit_length())
    span = min(BUCKET_BITS, 64 - 2 * bits)
    if span < 1:
        raise ValueError('Too many nodes to pack a face into 64 bits')
    # Sorted node ids of face k of tet t live at flat position 4t + k.
    lo, mid, hi = (np.empty((count, 4), dtype=np.uint32 if bits <= 32 else np.uint64)
                   for _ in range(3))
    for k, (a, b, c) in enumerate(FACE_NODES):
        x, y, z = tets[:, a], tets[:, b], tets[:, c]
        small = np.minimum(np.minimum(x, y), z)
        large = np.maximum(np.maximum(x, y), z)
        lo[:, k], hi[:, k] = small, large
        mid[:, k] = x + y + z - small - large
    lo, mid, hi = lo.reshape(-1), mid.reshape(-1), hi.reshape(-1)
    bucket = (lo >> span).astype(np.uint16 if bits - span <= 16 else np.uint32)
    order = np.argsort(bucket, kind='stable')
    cuts = np.searchsorted(bucket[order], np.arange(int(bucket.max()) + 2))
    del bucket
    flat = nbr.reshape(-1)

    def match(i):
        faces = order[cuts[i]:cuts[i + 1]]
        key = (((lo[faces] - (i << span)).astype(np.uint64) << np.uint64(2 * bits))
               | (mid[faces].astype(np.uint64) << np.uint64(bits)) | hi[faces])
        ix = np.argsort(key)
        faces, key = faces[ix], key[ix]
        same = key[1:] == key[:-1]
        if np.any(same[:-1] & same[1:]):
            raise ValueError('Non-manifold volume mesh: more than two tets share a face')
        pair = np.flatnonzero(same)
        a, b = faces[pair], faces[pair + 1]
        flat[a] = b // 4
        flat[b] = a // 4

    with ThreadPoolExecutor(max(1, int(threads))) as pool:
        list(pool.map(match, range(len(cuts) - 1)))
    return nbr


def subset_neighbours(nbr, keep):
    """Adjacency of tets[keep]; a face shared with a dropped tet becomes boundary."""
    index = np.full(len(nbr) + 1, -1, dtype=nbr.dtype)  # index[-1] maps -1 to -1
    index[np.flatnonzero(keep)] = np.arange(int(np.count_nonzero(keep)), dtype=nbr.dtype)
    return index[nbr[keep]]


def face_chunks(tets, nbr, chunk=2_000_000):
    """Yield (faces (m, 3), owner, other) covering every face exactly once.

    other is -1 for a boundary face; an interior face is reported by the
    lower-numbered of its two tets. Face nodes keep the tet's own ordering.
    """
    for start in range(0, len(tets), chunk):
        part = nbr[start:start + chunk]
        own = np.arange(start, start + len(part))[:, None]
        t, k = np.nonzero((part < 0) | (part > own))
        other = part[t, k].astype(np.int64)
        t += start
        yield tets[t[:, None], _FACE_NODES[k]], t, other
