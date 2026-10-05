"""Precompute exact vector cross-sections for panel A of the poster figure.

Replaces the old rasterised point-in-shell mask. Each tissue body and each lead
part is intersected with the plane z = z0 and chained into closed polygons; the
result is cached so the figure script stays fast and self-contained.

Validated against a 0.05 mm point-in-shell raster: all ten tissue/level
combinations agree to 0.000% by area (even-odd rule).

The plane is taken at the centre of the DRIVEN SOURCE contacts (local contact 4,
both rails), not at the midpoint of the whole driven span -- that midpoint falls
in the 1 mm gap between contacts 4 and 5, where no electrode exists to draw.

Run:  fem/.venv/bin/python src/neuron/build_cross_sections.py
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))
sys.path.insert(0, HERE)

import config as C                      # noqa: E402
from stlio import read_stl              # noqa: E402
from cross_section import section       # noqa: E402

OUT = os.path.join(HERE, "out", "poster_matched", "sections.npz")
TISSUES = ["epidural", "dura", "csf", "white", "grey"]
RUNS = {"dorsal": "63c73a9ab23d01d11669c7b4", "ventral": "1b4c9eec2c617fe8e4404f4a"}
# Centre of driven source contacts 4/12, measured from each run's lead STLs.
SLICE_Z = {"dorsal": 78.45, "ventral": 136.85}
DRIVEN = ("4", "12")                    # source rail; 5/13 are the sink, 4 mm rostral


def main():
    store, index = {}, {"slice_z": SLICE_Z, "tissues": TISSUES, "loops": {}}
    for model, z in SLICE_Z.items():
        for t in TISSUES:
            n = 0
            for f in C.TISSUE_BODIES[t]:
                v, tri = read_stl(f)
                for L in section(v, tri, z):
                    store["%s|%s|%d" % (model, t, n)] = L
                    n += 1
            index["loops"]["%s|%s" % (model, t)] = n
            print("  %-8s %-9s %2d loops" % (model, t, n), flush=True)

        lead = os.path.join(ROOT, "fem", "out", "lead_runs", RUNS[model], "lead")
        for role, files in (
                ("contact", [p for p in glob.glob(lead + "/SCS Lead Electrode *.stl")
                             if os.path.basename(p).split()[-1][:-4] in DRIVEN]),
                ("insulator", sorted(glob.glob(lead + "/insulator_lead_*.stl")))):
            n = 0
            for f in sorted(files):
                v, tri = read_stl(f)
                for L in section(v, tri, z):
                    store["%s|lead_%s|%d" % (model, role, n)] = L
                    n += 1
            index["loops"]["%s|lead_%s" % (model, role)] = n
            print("  %-8s lead_%-4s %2d loops" % (model, role, n), flush=True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez_compressed(OUT, index_json=json.dumps(index), **store)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
