"""Verify the P1 current-conduction solver against an exact analytic solution.

Problem: a conducting sphere shell a <= r <= b, conductivity sigma, current I
injected uniformly over r = a and the outer surface r = b held at 0 V.  Exactly

    V(r) = I / (4 pi sigma) * (1/r - 1/b)

This exercises the whole chain the real solve uses -- element assembly, the
mm -> m unit conversion, volumetric current injection, the Dirichlet pin -- so
if the units were wrong by 1000x it would show here and not in a plausible
looking picture of the cord.

Run:  python fem/scripts/verify_solver.py
"""
import os
import sys

import numpy as np
import scipy.sparse.linalg as spl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assign_and_solve import assemble, boundary_nodes

A_MM, B_MM, SIGMA, I = 4.0, 24.0, 0.3, 1.0


def main(h=0.30, grow=0.10):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("shell")
    outer = gmsh.model.occ.addSphere(0, 0, 0, B_MM)
    inner = gmsh.model.occ.addSphere(0, 0, 0, A_MM)
    gmsh.model.occ.cut([(3, outer)], [(3, inner)])
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    f = gmsh.model.mesh.field
    f.add("MathEval", 1)
    f.setString(1, "F", "%f + %f*(sqrt(x*x+y*y+z*z) - %f)" % (h, grow, A_MM))
    f.setAsBackgroundMesh(1)
    gmsh.model.mesh.generate(3)
    nt, nc, _ = gmsh.model.mesh.getNodes()
    o = np.argsort(nt)
    nodes = nc.reshape(-1, 3)[o]
    remap = np.zeros(int(nt.max()) + 1, dtype=np.int64)
    remap[nt[o]] = np.arange(len(nt))
    ety, _, en = gmsh.model.mesh.getElements(3)
    tets = remap[np.concatenate([e for t, e in zip(ety, en) if t == 4])].reshape(-1, 4)
    gmsh.finalize()

    sig = np.full(len(tets), SIGMA)
    K, vol, tets = assemble(nodes, tets, sig)
    r = np.linalg.norm(nodes, axis=1)

    bf = boundary_nodes(tets)
    bn = np.unique(bf)
    inner_n = bn[r[bn] < (A_MM + B_MM) / 2]
    outer_n = bn[r[bn] >= (A_MM + B_MM) / 2]

    # inject I over the inner surface, area-weighted
    p = nodes[bf] * 1e-3
    ar = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
    w = np.zeros(len(nodes))
    np.add.at(w, bf.ravel(), np.repeat(ar / 3.0, 3))
    src = np.zeros(len(nodes))
    src[inner_n] = w[inner_n] / w[inner_n].sum() * I

    free = np.ones(len(nodes), bool)
    free[outer_n] = False
    fi = np.flatnonzero(free)
    V = np.zeros(len(nodes))
    V[fi] = spl.spsolve(K[fi][:, fi].tocsc(), src[fi])

    exact = I / (4 * np.pi * SIGMA) * (1.0 / (r * 1e-3) - 1.0 / (B_MM * 1e-3))
    m = (r > A_MM * 1.25) & (r < B_MM * 0.95)
    err = np.abs(V[m] - exact[m])
    scale = exact[m].max()
    print("nodes %d tets %d" % (len(nodes), len(tets)))
    print("V at r=a : FEM %.4f V   exact %.4f V" % (V[inner_n].mean(), exact[inner_n].mean()))
    print("interior rel L-inf error %.3f %%   rel RMS %.3f %%"
          % (100 * err.max() / scale, 100 * np.sqrt((err ** 2).mean()) / scale))
    net = K.dot(V)
    print("net current in  %+.6f A   out %+.6f A"
          % (net[inner_n].sum(), net[outer_n].sum()))
    ok = err.max() / scale < 0.03 and abs(net[outer_n].sum() + I) < 1e-6
    print("VERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    codes = []
    for h, grow in [(0.60, 0.16), (0.40, 0.11), (0.28, 0.08)]:
        print("\n=== target h at r=a: %.2f mm, growth %.2f ===" % (h, grow))
        codes.append(main(h, grow))
    sys.exit(min(codes))
