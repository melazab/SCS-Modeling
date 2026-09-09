#!/usr/bin/env python3
"""Find electrically isolated bodies in the RADO-SCS STL set.

The full-model solve fails with

    Distributed sparse solver minimum pivot in absolute value =
    1.942856684E-07 at node 8711510 VOLT.
    *** ERROR *** There is at least 1 small equation solver pivot term ...
    Please check for an insufficiently constrained model.

which is the signature of a region with no conduction path to ground. The 245
bodies are joined only by bonded electric contact within a pinball radius, so
any body further than that from every neighbour is an island.

This works it out from the geometry directly: build a proximity graph over the
STL surfaces, find its connected components, and report any component that is
not the main one. Pure numpy -- scipy is not installed here.

    python3 find_floating_bodies.py --tol 0.3        # tolerance in mm

Tolerance is the surface-to-surface gap below which two bodies are taken to be
in contact. Sweep it (--sweep) to see how the island count depends on it, which
is exactly the sensitivity that the contact pinball radius controls in Ansys.
"""
import argparse
import os
import struct
import sys

import numpy as np


def read_stl_vertices(path):
    """Return an (N,3) float32 array of triangle vertices from a binary or ASCII STL."""
    with open(path, "rb") as fh:
        head = fh.read(84)
        if len(head) < 84:
            return np.empty((0, 3), np.float32)
        ntri = struct.unpack("<I", head[80:84])[0]
        rest = fh.read()
    if len(rest) == ntri * 50 and ntri > 0:            # binary
        rec = np.frombuffer(rest, dtype=np.uint8).reshape(ntri, 50)
        verts = rec[:, 12:48].copy().view("<f4").reshape(ntri * 3, 3)
        return np.asarray(verts, dtype=np.float32)
    # ASCII fallback
    pts = []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if "vertex" in line:
                pts.append([float(v) for v in line.split()[1:4]])
    return np.asarray(pts, dtype=np.float32)


def subsample(v, n):
    if len(v) <= n:
        return v
    idx = np.linspace(0, len(v) - 1, n).astype(np.int64)
    return v[idx]


def min_gap(a, b, chunk=512):
    """Smallest distance between two point sets, in chunks to bound memory."""
    best = np.inf
    for i in range(0, len(a), chunk):
        blk = a[i:i + chunk]
        d2 = ((blk[:, None, :] - b[None, :, :]) ** 2).sum(-1).min()
        if d2 < best:
            best = d2
            if best == 0.0:
                return 0.0
    return float(np.sqrt(best))


class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def components(names, boxes, pts, tol):
    n = len(names)
    dsu = DSU(n)
    lo = np.array([b[0] for b in boxes])
    hi = np.array([b[1] for b in boxes])
    pairs = 0
    for i in range(n):
        # bounding boxes must overlap once inflated by tol, else skip the pair
        cand = np.where((lo[i + 1:] <= hi[i] + tol).all(1) &
                        (hi[i + 1:] >= lo[i] - tol).all(1))[0] + i + 1
        for j in cand:
            if dsu.find(i) == dsu.find(j):
                continue
            pairs += 1
            if min_gap(pts[i], pts[j]) <= tol:
                dsu.union(i, j)
    groups = {}
    for i in range(n):
        groups.setdefault(dsu.find(i), []).append(i)
    return sorted(groups.values(), key=len, reverse=True), pairs


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--stl-dir", default=os.path.join(here, "..", "STL_files"))
    ap.add_argument("--tol", type=float, default=0.3, help="contact gap, mm")
    ap.add_argument("--points", type=int, default=1500, help="vertices sampled per body")
    ap.add_argument("--sweep", action="store_true", help="try a range of tolerances")
    args = ap.parse_args()

    files = sorted(f for f in os.listdir(args.stl_dir) if f.lower().endswith(".stl"))
    names, boxes, pts = [], [], []
    for f in files:
        v = read_stl_vertices(os.path.join(args.stl_dir, f))
        if len(v) == 0:
            print("  ! no vertices: %s" % f)
            continue
        names.append(os.path.splitext(f)[0])
        boxes.append((v.min(0), v.max(0)))
        pts.append(subsample(np.unique(v, axis=0), args.points))
    print("loaded %d bodies, %d sampled points total\n"
          % (len(names), sum(len(p) for p in pts)))

    tols = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0] if args.sweep else [args.tol]
    for tol in tols:
        comps, pairs = components(names, boxes, pts, tol)
        iso = [c for c in comps[1:]]
        print("tol %.2f mm -> %d component(s); main has %d bodies, %d isolated in %d group(s) "
              "(%d pairs tested)"
              % (tol, len(comps), len(comps[0]), sum(len(c) for c in iso), len(iso), pairs))
        if not args.sweep:
            print()
            for k, c in enumerate(comps):
                tag = "MAIN" if k == 0 else "ISLAND %d" % k
                print("%s -- %d bodies" % (tag, len(c)))
                if k > 0 or len(c) <= 12:
                    for i in sorted(c, key=lambda i: names[i]):
                        print("    %s" % names[i])
            if len(comps) > 1:
                print("\nThese islands have no conduction path to the rest of the model.")
                print("In Ansys each one contributes a near-zero pivot, which is the")
                print("'insufficiently constrained model' error. Options: ground them,")
                print("widen the contact pinball radius, or add the missing soft-tissue")
                print("envelope so they have something to conduct into.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
