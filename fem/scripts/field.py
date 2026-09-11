"""Sample the FEM potential at arbitrary points, and export a structured grid.

Interface note for the NEURON stage (docs/neuron_plan.md chunk 3): this module
deliberately matches that contract's shape -- origin / spacing / shape / one
`phi` per contact / explicit units -- WITHOUT implementing its unit convention.
Everything here is in the model's own units:

    coordinates   mm   (document coordinates of NBF_RADO-SCS; identity placements)
    phi           V/A  (transfer impedance, ohms) -- multiply by contact currents
    V             V    (the bipolar 1 A case, = phi[2] - phi[4])

The NEURON side converts once, in its own field.py: mm -> um, V -> mV.
"""
import numpy as np
from scipy.spatial import cKDTree


class TetField:
    def __init__(self, nodes, tets, values):
        """values: (nfield, nnode) nodal scalars."""
        self.nodes = np.ascontiguousarray(nodes, dtype=np.float64)
        self.tets = np.ascontiguousarray(tets, dtype=np.int64)
        self.values = np.atleast_2d(values)
        self.cen = self.nodes[self.tets].mean(axis=1)
        self.tree = cKDTree(self.cen)

    def locate(self, pts, kmax=96):
        """Return (tet index, barycentric weights); tet = -1 if outside."""
        pts = np.atleast_2d(np.asarray(pts, dtype=np.float64))
        out = np.full(len(pts), -1, dtype=np.int64)
        bw = np.zeros((len(pts), 4))
        todo = np.arange(len(pts))
        k = 8
        while todo.size and k <= kmax:
            _, cand = self.tree.query(pts[todo], k=k, workers=-1)
            cand = np.atleast_2d(cand)
            found = np.zeros(todo.size, dtype=bool)
            for c in range(cand.shape[1] if k == 8 else cand.shape[1]):
                if found.all():
                    break
                sel = ~found
                ti = cand[sel, c]
                w = self._bary(pts[todo[sel]], ti)
                ok = (w >= -1e-9).all(axis=1)
                gi = np.flatnonzero(sel)[ok]
                out[todo[gi]] = ti[ok]
                bw[todo[gi]] = w[ok]
                found[gi] = True
            todo = todo[~found]
            k *= 3
        return out, bw

    def _bary(self, p, ti):
        n = self.nodes[self.tets[ti]]
        T = np.stack([n[:, 1] - n[:, 0], n[:, 2] - n[:, 0], n[:, 3] - n[:, 0]], axis=2)
        det = np.linalg.det(T)
        det = np.where(np.abs(det) < 1e-20, 1e-20, det)
        rhs = (p - n[:, 0])[:, :, None]
        lam = np.linalg.solve(T, rhs)[:, :, 0]
        return np.concatenate([1 - lam.sum(axis=1, keepdims=True), lam], axis=1)

    def __call__(self, pts, fill=np.nan):
        ti, w = self.locate(pts)
        res = np.full((len(self.values), len(ti)), fill, dtype=np.float64)
        ok = ti >= 0
        if ok.any():
            vn = self.values[:, self.tets[ti[ok]]]          # (nf, nok, 4)
            res[:, ok] = np.einsum("fnk,nk->fn", vn, w[ok])
        return res


def load(path):
    d = np.load(path, allow_pickle=True)
    return d


def structured_export(sol_npz, out_npz, origin, spacing, shape):
    d = np.load(sol_npz, allow_pickle=True)
    f = TetField(d["nodes"], d["tets"], d["phi"].astype(np.float64))
    ax = [origin[i] + spacing * np.arange(shape[i]) for i in range(3)]
    gx, gy, gz = np.meshgrid(*ax, indexing="ij")
    P = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    phi = np.full((8, len(P)), np.nan)
    for s in range(0, len(P), 200_000):
        sl = slice(s, s + 200_000)
        phi[:, sl] = f(P[sl])
    inside = np.isfinite(phi[0])
    # tissue label on the same grid, for the NEURON stage to pick dorsal-column points
    ti, _ = f.locate(P)
    lab = np.where(ti >= 0, d["label"][np.clip(ti, 0, None)], -1).astype(np.int16)
    np.savez_compressed(
        out_npz,
        origin=np.asarray(origin, float), spacing=float(spacing),
        shape=np.asarray(shape, int),
        phi=phi.reshape((8,) + tuple(shape)).astype(np.float32),
        label=lab.reshape(shape), order=d["order"],
        units_coords="mm", units_phi="V/A (transfer impedance, ohm)",
        note=("phi[i] is the potential for +1 A into contact i+1 and -1 A spread "
              "over the insulating outer boundary; any zero-net-current montage "
              "is sum_i I_i*phi[i]. Bipolar 1 A src=3 sink=5 is phi[2]-phi[4]."),
        coverage=float(inside.mean()))
    return inside.mean()
