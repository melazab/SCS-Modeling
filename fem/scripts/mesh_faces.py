"""Match tetrahedron faces with bounded sort memory.

Large meshes are partitioned to disk by a hash of the sorted vertex IDs.
Shared faces always land in the same bucket. Scratch is on the output disk,
not /tmp (which is RAM-backed on this workstation).
"""
import os
import tempfile
import numpy as np
import config as C

FACE_NODES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def _match(rows):
    ix = np.lexsort((rows[:, 2], rows[:, 1], rows[:, 0]))
    rows = rows[ix]
    same = np.all(rows[1:, :3] == rows[:-1, :3], axis=1)
    if np.any(same[:-1] & same[1:]):
        raise ValueError('Non-manifold volume mesh: more than two tets share a face')
    paired = np.flatnonzero(same)
    lone = np.ones(len(rows), bool)
    lone[paired] = False
    lone[paired + 1] = False
    return rows[:, :3], rows[:, 3], paired, lone


def face_batches(tets, target_faces=1_000_000, chunk_tets=100_000):
    if len(tets) == 0:
        return
    dtype = np.uint32 if max(int(tets.max()), len(tets)) < 2**32 else np.int64
    buckets = max(1, int(np.ceil(4 * len(tets) / target_faces)))

    def records(start, part):
        for face in FACE_NODES:
            rows = np.empty((len(part), 4), dtype=dtype)
            rows[:, :3] = np.sort(part[:, face], axis=1)
            rows[:, 3] = np.arange(start, start + len(part))
            yield rows

    if buckets == 1:
        yield _match(np.concatenate(list(records(0, tets))))
        return
    os.makedirs(C.OUT, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='face-sort-', dir=C.OUT) as scratch:
        paths = [os.path.join(scratch, str(i)) for i in range(buckets)]
        for start in range(0, len(tets), chunk_tets):
            for rows in records(start, tets[start:start + chunk_tets]):
                key = (rows[:, 0].astype(np.uint64) * 73856093
                       ^ rows[:, 1].astype(np.uint64) * 19349663
                       ^ rows[:, 2].astype(np.uint64) * 83492791) % buckets
                ix = np.argsort(key)
                rows, key = rows[ix], key[ix]
                cuts = np.r_[0, np.flatnonzero(key[1:] != key[:-1]) + 1, len(key)]
                for a, b in zip(cuts[:-1], cuts[1:]):
                    with open(paths[int(key[a])], 'ab') as f:
                        rows[a:b].tofile(f)
            if start % (10 * chunk_tets) == 0:
                print('  face partition %d/%d tets' % (start, len(tets)), flush=True)
        for i, path in enumerate(paths):
            if os.path.exists(path):
                rows = np.fromfile(path, dtype=dtype).reshape(-1, 4)
                yield _match(rows)
                os.unlink(path)
                print('  face sort %d/%d partitions' % (i + 1, buckets), flush=True)
