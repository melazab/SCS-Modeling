"""Solve per-contact voltage basis fields on a completed live-lead mesh.

CLI: solve_lead.py <run_dir> --lead-dir <run_dir>/lead [--background]
Validates mesh/lead provenance, checks the memory budget, runs the shared P1
solver and atomically publishes an accepted solution and its manifest.
SCS Potential Visualizer owns this worker; SCS Mesh Generator creates its mesh.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from elmer_backend import solve as solve_elmer
from solve_resources import Resources


import artifacts


def main():
    # Own process group lets Abort stop MPI launchers and ranks as well.
    if hasattr(os, 'setsid') and os.getpgrp() != os.getpid():
        os.setsid()
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("lead_run_dir",
                   help="scratch directory already containing mesh.npz "
                        "(Stage 1's own output dir, e.g. "
                        "fem/out/lead_runs/<lead_name>/)")
    p.add_argument("--lead-dir", default=None,
                   help="lead export directory whose contact/insulator STLs "
                        "this mesh was built against; omit for the frozen "
                        "fem/leads/fem_dorsal_T10 lead")
    p.add_argument("--background", action="store_true",
                   help="keep the bone-conductivity background tets (matches "
                        "SCS_BACKGROUND=1 for assign_and_solve.py's main())")
    p.add_argument('--cpus', type=int, default=1, help='Local MPI tasks, one thread per task')
    p.add_argument('--memory-gb', type=float, default=0., help='Worker process-tree RSS budget in GB; 0 = automatic')
    args = p.parse_args()

    lead_run_dir = os.path.abspath(args.lead_run_dir)
    if lead_run_dir == os.path.abspath(C.OUT):
        raise SystemExit("refusing to solve into fem/out/ itself -- that is "
                          "the field of record; use a fem/out/lead_runs/<name>/ "
                          "scratch directory instead")

    mesh_npz = os.path.join(lead_run_dir, "mesh.npz")
    if not os.path.isfile(mesh_npz):
        raise SystemExit("no mesh.npz in %s -- run Stage 1 (Generate Mesh) "
                          "for this lead first" % lead_run_dir)

    if args.lead_dir is not None and not C.is_single_lead_dir(args.lead_dir):
        raise SystemExit(
            "--lead-dir %r is not a single-lead flat-layout export (missing "
            "its insulator STL, or has zero 'SCS Lead Electrode *.stl' "
            "contact files) -- a multi-lead export directory is out of "
            "scope here" % args.lead_dir)
    contact_stl, insulator_stl = C.lead_stls(args.lead_dir)

    if not artifacts.preview_ready(lead_run_dir):
        raise SystemExit('Mesh preview is incomplete or stale; run Generate Mesh first')
    with open(os.path.join(lead_run_dir, 'mesh_report.json')) as fh:
        mesh_report = json.load(fh)
    expected = artifacts.mesh_signature(mesh_report['params'], contact_stl, insulator_stl)
    if not artifacts.valid(mesh_npz, expected):
        raise SystemExit('Lead geometry does not match this mesh; refusing to solve')

    signature = artifacts.solution_signature(lead_run_dir, args.background)
    t0 = time.time()
    def progress(pct):
        print("PROGRESS %d" % pct, flush=True)
    with Resources(args.cpus, args.memory_gb) as resources:
        with np.load(mesh_npz) as d:
            nodes, tets = d['nodes'], d['tets']
        print('Elmer local solve: %d CPUs, %.2f GB budget; %d nodes, %d tets' %
              (resources.cpus, resources.memory_gb, len(nodes), len(tets)), flush=True)
        result = solve_elmer(nodes, tets, contact_stl, insulator_stl,
                             args.background, os.path.join(lead_run_dir, 'elmer_bg' if args.background else 'elmer'),
                             resources, progress)
        resource_records = resources.records

        out_name = "solution_bg.npz" if args.background else "solution.npz"
        out_path = os.path.join(lead_run_dir, out_name)
        if signature != artifacts.solution_signature(lead_run_dir, args.background):
            raise RuntimeError("Model inputs changed during solve; result not published")
        if expected != artifacts.mesh_signature(mesh_report['params'], contact_stl, insulator_stl):
            raise RuntimeError('Lead inputs changed during solve; result not published')
        artifacts.atomic_npz(out_path, **result)
        artifacts.publish(out_path, signature)
    elapsed = time.time() - t0
    print("wrote %s  %.0fs" % (out_name, elapsed))

    phi = result["phi"]
    report = dict(
        lead_run_dir=lead_run_dir,
        lead_dir=args.lead_dir,
        solution_npz=out_path,
        nodes=int(len(result["nodes"])), tets=int(len(result["tets"])),
        background_included=bool(args.background),
        elapsed_s=round(elapsed, 1),
        backend=str(result['backend']), cpus=int(result['cpus']),
        memory_budget_gb=float(result['memory_budget_gb']),
        peak_rss_gb=float(result['peak_rss_gb']), resource_stages=resource_records,
        relative_residuals=result["relative_residuals"].tolist(),
        accepted_residual=float(result["accepted_residual"]),
        phi_min=float(phi.min()), phi_max=float(phi.max()),
        v_min=float(result["V"].min()), v_max=float(result["V"].max()),
    )
    print("SOLVE_JSON " + json.dumps(report))


if __name__ == "__main__":
    main()
