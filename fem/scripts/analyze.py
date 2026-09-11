"""Per-tissue field statistics, and the one comparison the paper lets us make.

Khadka et al. 2020 report, for their exemplary bipolar 1 A simulation:
    peak electric field   12 kV/m in white matter, 4.2 kV/m in grey matter
    peak surface voltage  1.2 kV
Those are the only published numbers this stripped model can be held against.
They come from the FULL model (19 compartments, ~150M elements, anisotropic
white matter), so agreement to a factor of ~2 would already be encouraging and
exact agreement would be suspicious.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C


def elem_gradients(nodes, tets):
    p = nodes[tets] * 1e-3
    J = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=2)
    vol = np.linalg.det(J) / 6.0
    Jinv = np.linalg.inv(J)
    g = np.empty((len(tets), 4, 3))
    g[:, 1:, :] = Jinv
    g[:, 0, :] = -g[:, 1:, :].sum(axis=1)
    return g, np.abs(vol)



def speckle_report(nodes, tets, lab, order):
    """How noisy the per-tet classification is.

    inside.py decides with a single +Z ray, and src/freecad/build_lead_config.py
    documents that single-ray parity is brittle: a ray grazing a shared edge can
    miscount one crossing and invert the answer. That agent votes across eight
    directions because it tests only 4284 lead-surface points; voting over
    millions of tet centroids would be four to eight times the cost of the whole
    classification. So instead of voting, this MEASURES the resulting noise:
    an isolated tetrahedron whose tissue differs from every one of its face
    neighbours is almost certainly a miscount, not anatomy. Isolated elements
    barely perturb a diffusion problem, but the rate is worth knowing and worth
    reporting rather than discovering later.
    """
    f = np.concatenate([tets[:, [0, 1, 2]], tets[:, [0, 1, 3]],
                        tets[:, [0, 2, 3]], tets[:, [1, 2, 3]]])
    owner = np.tile(np.arange(len(tets)), 4)
    key = np.sort(f, axis=1)
    o = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    key, owner = key[o], owner[o]
    same = np.all(key[1:] == key[:-1], axis=1)
    i0 = np.flatnonzero(same)
    a, b = owner[i0], owner[i0 + 1]
    nb = np.zeros(len(tets), dtype=np.int32)
    agree = np.zeros(len(tets), dtype=np.int32)
    for x, y in ((a, b), (b, a)):
        np.add.at(nb, x, 1)
        np.add.at(agree, x, (lab[x] == lab[y]).astype(np.int32))
    iso = (nb > 0) & (agree == 0)
    print("\nclassification noise")
    print("  isolated tets (tissue unlike every face neighbour): %d of %d = %.4f %%"
          % (iso.sum(), len(tets), 100 * iso.mean()))
    for i, name in enumerate(order):
        n = int((iso & (lab == i)).sum())
        if n:
            print("      %-11s %6d  (%.4f %% of that tissue)"
                  % (name, n, 100 * n / max((lab == i).sum(), 1)))

def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "solution.npz"
    print("=== %s ===" % name)
    d = np.load(os.path.join(C.OUT, name), allow_pickle=True)
    nodes, tets, lab = d["nodes"], d["tets"], d["label"]
    V = d["V"].astype(np.float64)
    order = [str(x) for x in d["order"]]
    g, vol = elem_gradients(nodes, tets)
    E = -np.einsum("eki,ek->ei", g, V[tets])       # V/m, constant per element
    Em = np.linalg.norm(E, axis=1)
    J = None

    print("bipolar %.0f A, +contact %d / -contact %d\n"
          % (C.CURRENT_A, C.SOURCE_CONTACT, C.SINK_CONTACT))
    print("%-11s %9s %11s %11s %11s %11s" %
          ("tissue", "tets", "vol mm3", "|E| max", "|E| p99.9", "|E| p99"))
    print("%-11s %9s %11s %11s %11s %11s" % ("", "", "", "kV/m", "kV/m", "kV/m"))
    for i, name in enumerate(order):
        m = lab == i
        if not m.any():
            continue
        e = Em[m] / 1e3
        print("%-11s %9d %11.2f %11.3f %11.3f %11.3f"
              % (name, m.sum(), vol[m].sum() * 1e9, e.max(),
                 np.percentile(e, 99.9), np.percentile(e, 99)))
    m = lab == -1
    if m.any():
        e = Em[m] / 1e3
        print("%-11s %9d %11.2f %11.3f %11.3f %11.3f"
              % ("background", m.sum(), vol[m].sum() * 1e9, e.max(),
                 np.percentile(e, 99.9), np.percentile(e, 99)))

    print("\nvoltage")
    print("  V max %+.1f V   V min %+.1f V   range %.1f V" % (V.max(), V.min(), V.max() - V.min()))
    for cid in (C.SOURCE_CONTACT, C.SINK_CONTACT):
        nd = np.unique(tets[lab == order.index("contact%d" % cid)])
        print("  contact %d mean %+.1f V" % (cid, V[nd].mean()))
    n3 = np.unique(tets[lab == order.index("contact%d" % C.SOURCE_CONTACT)])
    n5 = np.unique(tets[lab == order.index("contact%d" % C.SINK_CONTACT)])
    R = (V[n3].mean() - V[n5].mean()) / C.CURRENT_A
    print("  bipolar transfer impedance  %.1f ohm  (contact 3 -> contact 5)" % R)

    print("\nagainst Khadka 2020 (full 19-compartment model, bipolar 1 A):")
    for name, paper in (("white", 12.0), ("grey", 4.2)):
        i = order.index(name)
        mine = np.percentile(Em[lab == i] / 1e3, 99.9)
        print("  peak |E| %-6s  paper %5.1f kV/m   here %6.2f kV/m (p99.9)  ratio %.2f"
              % (name, paper, mine, mine / paper))
    print("  peak surface voltage  paper 1.2 kV   here %.2f kV (max |V|)"
          % (max(abs(V.max()), abs(V.min())) / 1e3))
    speckle_report(nodes, tets, lab, order)


if __name__ == "__main__":
    main()
