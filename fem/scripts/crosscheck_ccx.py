"""Independent cross-check of the P1 solver using CalculiX, via the heat analogy.

THE ANALOGY, STATED ONCE AND LOUDLY
-----------------------------------
Steady heat conduction and steady current conduction are the SAME equation.
CalculiX has no electric-conduction step, so the deck below is a heat-transfer
deck whose symbols mean electrical quantities:

    CalculiX symbol            what it means here          unit
    -------------------------  --------------------------  ------------
    NT  (nodal temperature)    electric potential V        volts
    *CONDUCTIVITY  k           electrical conductivity     S/m
    *CFLUX  (concentrated q)   injected current            amperes
    *BOUNDARY dof 11           prescribed potential        volts

Node coordinates are written in METRES so that k in S/m and flux in A give a
result in volts with no further scaling.  Anyone reading fem/out/ccx/*.frd must
read "temperature" as "volts".  This file exists ONLY as a second opinion on
fem/scripts/assign_and_solve.py -- the field of record is solution.npz, which is
in volts throughout and needs no translation.

Run after assign_and_solve.py:  python fem/scripts/crosscheck_ccx.py
"""
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

CCX = "/bin/ccx"


def write_deck(d, path):
    nodes, tets, lab = d["nodes"], d["tets"], d["label"]
    order = [str(x) for x in d["order"]]
    with open(path, "w") as f:
        f.write("** see fem/scripts/crosscheck_ccx.py: NT means VOLTS\n*NODE, NSET=Nall\n")
        np.savetxt(f, np.column_stack([np.arange(1, len(nodes) + 1), nodes * 1e-3]),
                   fmt="%d, %.9g, %.9g, %.9g")
        for u in sorted(set(int(x) for x in np.unique(lab))):
            name = order[u] if u >= 0 else "background"
            m = np.flatnonzero(lab == u)
            f.write("*ELEMENT, TYPE=C3D4, ELSET=E%s\n" % name)
            np.savetxt(f, np.column_stack([m + 1, tets[m] + 1]), fmt="%d, %d, %d, %d, %d")
            sig = (C.SIGMA_METAL if name.startswith("contact")
                   else C.SIGMA.get(name, C.SIGMA["background"]))
            f.write("*MATERIAL, NAME=M%s\n*CONDUCTIVITY\n%.8g\n" % (name, sig))
            f.write("*SOLID SECTION, ELSET=E%s, MATERIAL=M%s\n" % (name, name))
        f.write("*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNall, 0.\n")
        f.write("*STEP\n*HEAT TRANSFER, STEADY STATE, SOLVER=ITERATIVE CHOLESKY\n")
        f.write("*BOUNDARY\n%d, 11, 11, 0.\n" % (int(d["pinned_node"]) + 1))
        f.write("*CFLUX\n")
        for cid, amps in ((C.SOURCE_CONTACT, +C.CURRENT_A), (C.SINK_CONTACT, -C.CURRENT_A)):
            m = lab == order.index("contact%d" % cid)
            nd, w = np.unique(tets[m], return_counts=True)
            w = w / w.sum() * amps
            for k, v in zip(nd, w):
                f.write("%d, 11, %.10g\n" % (k + 1, v))
        f.write("*NODE FILE\nNT\n*END STEP\n")


def read_frd(path, nnode):
    V = np.full(nnode, np.nan)
    with open(path) as f:
        active = False
        for line in f:
            if line.startswith(" -4") and "NDTEMP" in line:
                active = True
                continue
            if active:
                if line.startswith(" -3"):
                    break
                if line.startswith(" -1"):
                    V[int(line[3:13]) - 1] = float(line[13:25])
    return V


def main():
    if not os.path.exists(CCX):
        raise SystemExit("ccx not found at %s" % CCX)
    d = np.load(os.path.join(C.OUT, "solution.npz"), allow_pickle=True)
    wd = os.path.join(C.OUT, "ccx")
    os.makedirs(wd, exist_ok=True)
    inp = os.path.join(wd, "scs")
    write_deck(d, inp + ".inp")
    print("deck written (%.0f MB); running ccx ..."
          % (os.path.getsize(inp + ".inp") / 1e6), flush=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(os.cpu_count()))
    r = subprocess.run([CCX, "scs"], cwd=wd, env=env, capture_output=True, text=True)
    tail = (r.stdout or "")[-1500:]
    print(tail)
    if not os.path.exists(inp + ".frd"):
        raise SystemExit("ccx produced no .frd\n" + (r.stderr or "")[-2000:])
    Vc = read_frd(inp + ".frd", len(d["nodes"]))
    Vp = d["V"].astype(np.float64)
    ok = np.isfinite(Vc)
    # both are defined up to a constant that the same pin removes, so compare directly
    rng = Vp.max() - Vp.min()
    diff = np.abs(Vc[ok] - Vp[ok])
    print("\nCROSS-CHECK  P1 (this repo) vs CalculiX heat analogy, same mesh")
    print("  nodes compared      %d" % ok.sum())
    print("  max |dV|            %.4e V  (%.4f %% of the %.3f V range)"
          % (diff.max(), 100 * diff.max() / rng, rng))
    print("  rms |dV|            %.4e V  (%.5f %%)"
          % (np.sqrt((diff ** 2).mean()), 100 * np.sqrt((diff ** 2).mean()) / rng))
    print("  VERDICT:", "AGREE" if diff.max() / rng < 0.01 else "DISAGREE")


if __name__ == "__main__":
    main()
