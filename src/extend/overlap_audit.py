"""Measure where anatomical bodies overlap each other.

CLI: fem/.venv/bin/python src/extend/overlap_audit.py [geometry_dir] [--h 0.25]
     [--workers N] [--json out.json]

Every closed anatomical STL in the directory (default: RADO's STL_files/) is
sampled on one shared grid of spacing h with the pipeline's own
point-in-shell test (fem/scripts/inside.py). A grid point inside two bodies of
DIFFERENT tissue classes is an overlap: the mesher can only give it one label,
so one of the two tissues silently loses that volume. Same-class overlaps are
reported separately; classification unions bodies of a class, so they cost
nothing in the FEM but are still double geometry.

Volumes are grid counts x h^3, so treat anything under a few h^3 as noise.
"""
import argparse
import json
import multiprocessing
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'fem', 'scripts'))
import config as C                                  # noqa: E402
from check_tissue_map import matches                # noqa: E402
from inside import InsideTester                     # noqa: E402
from stlio import read_stl                          # noqa: E402

_JOB = None


def classify(directory):
    """[(path, stem, pipeline class)] for every anatomical STL at the top level."""
    entries, _ = C._load_tissue_map()
    names = dict((cls, name) for name, cls in C.ANATOMY_CLASS)
    out = []
    for fn in sorted(os.listdir(directory)):
        if not fn.lower().endswith('.stl'):
            continue
        stem = os.path.splitext(fn)[0]
        hit = next((t['name'] for t in entries if any(matches(stem, p) for p in t['patterns'])), None)
        if hit in names:
            out.append((os.path.join(directory, fn), stem, names[hit]))
    return out


def _inside_indices(i):
    """Global grid indices inside body i (sorted int64)."""
    bodies, origin, h, shape = _JOB
    v, t = read_stl(bodies[i][0])
    lo = np.floor((v.min(0) - origin) / h).astype(int)
    hi = np.ceil((v.max(0) - origin) / h).astype(int)
    lo, hi = np.maximum(lo, 0), np.minimum(hi, np.array(shape) - 1)
    ax = [np.arange(lo[d], hi[d] + 1) for d in range(3)]
    tester = InsideTester(v, t)
    found = []
    for ix in ax[0]:   # one x-slab at a time keeps memory flat
        jy, kz = np.meshgrid(ax[1], ax[2], indexing='ij')
        pts = np.column_stack([np.full(jy.size, ix), jy.ravel(), kz.ravel()])
        hit = tester(origin + pts * h)
        if hit.any():
            p = pts[hit]
            found.append((p[:, 0] * shape[1] + p[:, 1]) * shape[2] + p[:, 2])
    return np.sort(np.concatenate(found)) if found else np.zeros(0, np.int64)


def audit(directory, h=0.25, workers=None):
    global _JOB
    bodies = classify(directory)
    mins, maxs = [], []
    for path, _, _ in bodies:
        v, _ = read_stl(path)
        mins.append(v.min(0)); maxs.append(v.max(0))
    mins, maxs = np.array(mins), np.array(maxs)
    origin = mins.min(0) - h
    shape = tuple(int(x) for x in np.ceil((maxs.max(0) + h - origin) / h) + 1)
    _JOB = (bodies, origin, h, shape)
    with multiprocessing.get_context('fork').Pool(workers or os.cpu_count()) as pool:
        idx = pool.map(_inside_indices, range(len(bodies)))
    _JOB = None
    voxel = h ** 3
    pairs = []
    for i in range(len(bodies)):
        for j in range(i + 1, len(bodies)):
            if np.any(maxs[i] < mins[j]) or np.any(maxs[j] < mins[i]):
                continue
            n = np.intersect1d(idx[i], idx[j], assume_unique=True).size
            if n:
                pairs.append(dict(a=bodies[i][1], b=bodies[j][1], class_a=bodies[i][2], class_b=bodies[j][2],
                                  overlap_mm3=n * voxel, frac_of_a=n / max(idx[i].size, 1),
                                  frac_of_b=n / max(idx[j].size, 1)))
    volumes = {bodies[i][1]: idx[i].size * voxel for i in range(len(bodies))}
    return bodies, volumes, pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('directory', nargs='?', default=C.STL)
    ap.add_argument('--h', type=float, default=0.25)
    ap.add_argument('--workers', type=int, default=None)
    ap.add_argument('--json', default=None)
    args = ap.parse_args()
    t0 = time.time()
    bodies, volumes, pairs = audit(args.directory, args.h, args.workers)
    cross = [p for p in pairs if p['class_a'] != p['class_b']]
    same = [p for p in pairs if p['class_a'] == p['class_b']]
    print('%d anatomical bodies in %s, grid %.2f mm, %.0f s' % (len(bodies), args.directory, args.h, time.time() - t0))
    print('different-class overlaps: %d pairs, %.1f mm3 in total' % (len(cross), sum(p['overlap_mm3'] for p in cross)))
    by = {}
    for p in cross:
        key = tuple(sorted((p['class_a'], p['class_b'])))
        n, v = by.get(key, (0, 0.0))
        by[key] = (n + 1, v + p['overlap_mm3'])
    for key, (n, v) in sorted(by.items(), key=lambda kv: -kv[1][1]):
        print('   %-12s x %-12s %4d pairs %9.1f mm3' % (key[0], key[1], n, v))
    print('largest:')
    for p in sorted(cross, key=lambda p: -p['overlap_mm3'])[:15]:
        print('   %8.1f mm3  %-45s (%s, %.1f%%)  x  %-45s (%s, %.1f%%)' % (
            p['overlap_mm3'], p['a'][8:][:45], p['class_a'], 100 * p['frac_of_a'],
            p['b'][8:][:45], p['class_b'], 100 * p['frac_of_b']))
    print('same-class overlaps: %d pairs, %.1f mm3' % (len(same), sum(p['overlap_mm3'] for p in same)))
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(dict(directory=args.directory, h=args.h, volumes_mm3=volumes, pairs=pairs), f, indent=1)


if __name__ == '__main__':
    main()
