"""Fast point-in-closed-triangle-mesh test (ray casting, XY-grid accelerated).

RADO's compartments are closed, manifold shells that tile space: each one is the
ACTUAL material region (the dura STL is the dura wall, not the sac interior), so
"inside body X" means "made of tissue X". Verified in fem/README.md.
"""
import numpy as np


class InsideTester:
    """Cast a +Z ray from each query point and count crossings (odd => inside)."""

    def __init__(self, verts, tris, ncell=96):
        self.v = np.ascontiguousarray(verts, dtype=np.float64)
        self.t = np.ascontiguousarray(tris, dtype=np.int64)
        self.a = self.v[self.t[:, 0]]
        self.b = self.v[self.t[:, 1]]
        self.c = self.v[self.t[:, 2]]
        self.lo = self.v.min(0)
        self.hi = self.v.max(0)
        span = np.maximum(self.hi[:2] - self.lo[:2], 1e-9)
        self.org = self.lo[:2] - 1e-6
        self.cell = span / ncell * (1 + 2e-6)
        self.n = ncell
        self._build()

    def _build(self):
        tlo = np.minimum(np.minimum(self.a, self.b), self.c)[:, :2]
        thi = np.maximum(np.maximum(self.a, self.b), self.c)[:, :2]
        i0 = np.clip(((tlo - self.org) / self.cell).astype(int), 0, self.n - 1)
        i1 = np.clip(((thi - self.org) / self.cell).astype(int), 0, self.n - 1)
        buckets = [[] for _ in range(self.n * self.n)]
        for k in range(len(self.t)):
            for ix in range(i0[k, 0], i1[k, 0] + 1):
                base = ix * self.n
                for iy in range(i0[k, 1], i1[k, 1] + 1):
                    buckets[base + iy].append(k)
        self.start = np.zeros(self.n * self.n + 1, dtype=np.int64)
        for i, bl in enumerate(buckets):
            self.start[i + 1] = self.start[i] + len(bl)
        self.items = np.fromiter((k for bl in buckets for k in bl),
                                 dtype=np.int64, count=int(self.start[-1]))

    def __call__(self, pts):
        pts = np.atleast_2d(np.asarray(pts, dtype=np.float64))
        out = np.zeros(len(pts), dtype=bool)
        inbox = np.all((pts >= self.lo - 1e-9) & (pts <= self.hi + 1e-9), axis=1)
        idx = np.flatnonzero(inbox)
        if idx.size == 0:
            return out
        ci = np.clip(((pts[idx, :2] - self.org) / self.cell).astype(int), 0, self.n - 1)
        flat = ci[:, 0] * self.n + ci[:, 1]
        order = np.argsort(flat, kind="stable")
        idx, flat = idx[order], flat[order]
        bnd = np.flatnonzero(np.diff(flat)) + 1
        for grp in np.split(np.arange(len(idx)), bnd):
            cellid = flat[grp[0]]
            tl = self.items[self.start[cellid]:self.start[cellid + 1]]
            if tl.size == 0:
                continue
            out[idx[grp]] = self._cross(pts[idx[grp]], tl)
        return out

    def _cross(self, p, tl):
        """Odd/even crossings of the +Z ray from p against triangles tl."""
        a, b, c = self.a[tl], self.b[tl], self.c[tl]
        # barycentric solve in XY for every (point, triangle) pair
        px = p[:, None, 0]
        py = p[:, None, 1]
        ax, ay = a[None, :, 0], a[None, :, 1]
        v0x, v0y = b[None, :, 0] - ax, b[None, :, 1] - ay
        v1x, v1y = c[None, :, 0] - ax, c[None, :, 1] - ay
        den = v0x * v1y - v1x * v0y
        safe = np.where(np.abs(den) < 1e-18, 1.0, den)
        qx, qy = px - ax, py - ay
        u = (qx * v1y - v1x * qy) / safe
        w = (v0x * qy - qx * v0y) / safe
        hit = (np.abs(den) >= 1e-18) & (u >= 0) & (w >= 0) & (u + w <= 1.0)
        zhit = (a[None, :, 2] + u * (b[None, :, 2] - a[None, :, 2])
                + w * (c[None, :, 2] - a[None, :, 2]))
        hit &= zhit > p[:, None, 2]
        return (hit.sum(axis=1) % 2) == 1


def classify(pts, testers, order):
    """First tester in `order` whose region contains the point wins.

    Returns an int array of indices into `order` (-1 = outside everything).
    """
    lab = np.full(len(pts), -1, dtype=np.int32)
    todo = np.arange(len(pts))
    for i, name in enumerate(order):
        if todo.size == 0:
            break
        m = testers[name](pts[todo])
        lab[todo[m]] = i
        todo = todo[~m]
    return lab
