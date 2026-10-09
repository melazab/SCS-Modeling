"""Keep one mesh_cost.CostModel alive and answer parameter sets over a pipe.

WHY A SERVER AND NOT A FUNCTION CALL

`fem/scripts/mesh_cost.py` exists so the SCS Mesh Generator panel can say
"3.5 M tets, ~519-864 MB RAM, ~3 min meshing" BEFORE Mohamed spends twenty
minutes finding out (see that module's docstring for the 2026-09-14 build
that implied ~226 M tets and could never have finished). Wiring it into the
panel runs into two hard constraints that together force this shape:

 1. FreeCAD's own Python has NO scipy. Every fem/scripts module the macros
    import keeps its scipy import inside a function for exactly this reason
    (build_mesh.build(), mesh_cost.CostModel.__init__). `import mesh_cost`
    works in FreeCAD; `mesh_cost.CostModel(...)` does not -- it needs
    scipy.spatial.cKDTree. So the model has to live in fem/.venv, i.e. in a
    subprocess.

 2. The expensive part is per-GEOMETRY, not per-parameter-set.
    CostModel.__init__ is ~5 s (reading the STLs, densifying them, and
    building one cKDTree per size-field region); CostModel.estimate() is
    ~1.5 ms. Re-running a one-shot script per spinbox change would pay the
    5 s every time and be useless as a live readout.

So: this process builds the model ONCE for one lead's geometry, prints a
`ready` line, and then answers one JSON request per stdin line with one JSON
response per stdout line, ~1.5 ms apart, until stdin closes or it is told to
QUIT. The panel keeps it alive for as long as the live lead's fingerprint is
unchanged and restarts it when the lead moves.

NOT AsyncProcessRunner. That class (fem/scripts/async_runner.py) runs a
one-shot script that prints PROGRESS lines and exits with a final JSON
report -- a completely different lifecycle from a long-lived request/response
pipe, and it has no way to write to a child's stdin. The panel drives this
one with a plain QProcess of its own.

PROTOCOL
    stdout, first line:  {"ready": true, "ram_gb": 30.4, "n_samples": 123456}
    stdin,  per request: {"id": 7, "H_MIN": 0.2, "H_MAX": 2.5,
                          "FIELDS": [["dura",0.25,0.4,1.0], ...]}
    stdout, per reply:   {"id": 7, "sev": "ok", "text": "3.5 M tets, ...",
                          "n_tets": 3500000.0, "ram_gb_lo": ..., "ram_gb_hi": ...,
                          "mesh_seconds": ...}
    stdin "QUIT" (or EOF) exits 0.
Any failure is reported as {"id": ..., "error": "..."} rather than a crash,
so a malformed request never takes the readout down with it.

CLI:
    python mesh_cost_server.py [<lead_dir>]
`lead_dir` is a single-lead flat-layout directory (config.is_single_lead_dir);
omit it, or pass "-", for the frozen fem/leads/fem_dorsal_T10 lead.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import mesh_cost


def _params_from(req):
    """build_mesh.DEFAULT_PARAMS-shaped dict out of one request. Only the
    four keys CostModel.raw_tets() actually reads are required; GRID/MARGIN
    are the build's own sampling/box settings and play no part in the
    estimate (the estimator has its own ESTIMATE_GRID_MM)."""
    return dict(H_MIN=float(req["H_MIN"]), H_MAX=float(req["H_MAX"]),
                FIELDS=[(str(f[0]), float(f[1]), float(f[2]), float(f[3]))
                        for f in req["FIELDS"]])


def main():
    lead_dir = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != "-" else None
    contact_stl, insulator_stl = C.lead_stls(lead_dir)
    model = mesh_cost.CostModel(contact_stl=contact_stl, insulator_stl=insulator_stl)
    ram_gb = mesh_cost.available_ram_gb()
    # `pid` is this process's OWN os.getpid(), reported because QProcess's
    # processId() was observed (live, 2026-09-15) to hand back a pid that
    # does not exist -- it can catch a transient fork during interpreter
    # startup rather than the process that ends up running this script.
    # CostEstimator.stop() needs a pid it can actually SIGKILL, so it gets
    # one straight from the horse's mouth.
    print(json.dumps(dict(ready=True, ram_gb=ram_gb, n_samples=model.n_samples,
                          pid=os.getpid())), flush=True)

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        if line == "QUIT":
            break
        rid = None
        try:
            req = json.loads(line)
            rid = req.get("id")
            est = model.estimate(_params_from(req))
            # A cluster job sends its allocation; otherwise this machine's RAM.
            sev, text = mesh_cost.describe(est, req.get("ram_gb"), req.get("budget", "available"))
            out = dict(id=rid, sev=sev, text=text)
            out.update({k: float(v) for k, v in est.items()})
        except Exception as exc:                                # noqa: BLE001
            out = dict(id=rid, error="%s: %s" % (type(exc).__name__, exc))
        print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
