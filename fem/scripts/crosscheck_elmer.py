"""Independent THIRD solve of the same problem, via Elmer's own StatCurrentSolver
-- no analogy this time (unlike crosscheck_ccx.py's heat-transfer stand-in),
this equation IS div(sigma grad V) = 0 in Elmer's own electric-conduction
solver, with no code shared with fem/scripts/assign_and_solve.py.

WHAT THIS SCRIPT DOES
----------------------
1. Reads fem/out/solution.npz for the two numbers Elmer needs that the P1
   solve already measured on this exact mesh: the volume of contact3 and
   contact5 (to convert +-1 A into a volumetric Current Source in A/m^3) and
   the pinned gauge node (so both solutions share the same V=0 reference and
   are directly comparable, not just comparable up to an additive constant).
2. Writes fem/out/elmer/case.sif and runs `ElmerGrid` (gmsh -> Elmer format,
   once) then `ElmerSolver`, both from the /opt/elmerfem install.
3. Parses the resulting raw-appended-binary VTU by hand (numpy only -- this
   box has no meshio/vtk in the system Python) and reports the same
   per-tissue |E| table and P1-vs-Elmer node-by-node voltage diff that
   crosscheck_ccx.py reports for CalculiX.

Boundary conditions are IDENTICAL to assign_and_solve.py's, restated in
Elmer's terms: contacts are sigma=1e4 S/m bodies (SIGMA_METAL, same clamp);
the driven pair gets a uniform volumetric Current Source over its own body
so the integral is exactly +-1 A; the six others get no source, which is an
exact floating conductor by conservation; the outer canal wall is Elmer's
default zero-flux Neumann BC (no keyword needed -- the mesh carries no 2D
physical surfaces to attach one to anyway); the all-Neumann gauge is fixed
by a single-node Dirichlet BC at the same node the P1 solve pins.

Run after export_results.py has written fem/out/mesh_tagged.msh:
    python fem/scripts/crosscheck_elmer.py
"""
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

ELMER_HOME = "/opt/elmerfem"
ELMER_BIN = os.path.join(ELMER_HOME, "bin")
ELMER_ENV = dict(
    os.environ,
    PATH=ELMER_BIN + os.pathsep + os.environ.get("PATH", ""),
    ELMER_HOME=ELMER_HOME,
    ELMER_SOLVER_HOME=os.path.join(ELMER_HOME, "share", "elmersolver"),
)

WD = os.path.join(C.OUT, "elmer")
MESH_MSH = os.path.join(C.OUT, "mesh_tagged.msh")
MESH_DIR = os.path.join(WD, "mesh_tagged")

SIF_TEMPLATE = """\
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "mesh_tagged"
  Include Path ""
  Results Directory "results"
End

Simulation
  Max Output Level = 5
  Coordinate System = Cartesian 3D
  Coordinate Scaling = 1.0e-3
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Post File = "case.vtu"
  Output File = "case.result"
End

Constants
  Permittivity of Vacuum = 8.8542e-12
End

{bodies}

{materials}

{body_forces}

Equation 1
  Name = "StatCurrent"
  Active Solvers(1) = 1
End

Solver 1
  Equation = Stat Current Solver
  Procedure = "StatCurrentSolve" "StatCurrentSolver"
  Variable = Potential
  Variable DOFs = 1
  Calculate Volume Current = True
  Calculate Joule Heating = False
  Linear System Solver = Direct
  Linear System Direct Method = MUMPS
  Steady State Convergence Tolerance = 1.0e-8
End

Boundary Condition 1
  Name = "gauge_pin"
  Target Nodes(1) = {pin_node}
  Potential = 0.0
End
"""

# ORDER, tags 1..len(C.ORDER) -- must match export_results.write_tagged_msh's
# numbering. That was 14 tags for the five-compartment model and is 20 now that
# config.ORDER carries all eleven anatomy classes; both ends derive it from
# C.ORDER, so nothing here had to change when the tissue set grew.
TAGS = list(range(1, len(C.ORDER) + 1))
MAT_OF = {}  # tissue name -> material index (metal contacts share one)
MATERIALS = []  # (name, sigma)
for name in C.ORDER:
    sigma = C.SIGMA_METAL if name.startswith("contact") else C.SIGMA[name]
    key = "metal_clamped" if name.startswith("contact") else name
    if key not in MAT_OF:
        MAT_OF[key] = len(MATERIALS) + 1
        MATERIALS.append((key, sigma))


def build_sif(vol3_m3, vol5_m3, pin_node_1idx):
    bodies = []
    for tag, name in zip(TAGS, C.ORDER):
        key = "metal_clamped" if name.startswith("contact") else name
        extra = ""
        if name == "contact%d" % C.SOURCE_CONTACT:
            extra = "\n  Body Force = 1"
        elif name == "contact%d" % C.SINK_CONTACT:
            extra = "\n  Body Force = 2"
        bodies.append('Body %d\n  Name = "%s"\n  Target Bodies(1) = %d\n'
                       '  Equation = 1\n  Material = %d%s\nEnd'
                       % (tag, name, tag, MAT_OF[key], extra))
    materials = []
    for i, (key, sigma) in enumerate(MATERIALS, start=1):
        materials.append('Material %d\n  Name = "%s"\n  Electric Conductivity = %.8g\nEnd'
                          % (i, key, sigma))
    j3 = C.CURRENT_A / vol3_m3
    j5 = -C.CURRENT_A / vol5_m3
    body_forces = (
        'Body Force 1\n  Name = "inject_source"\n  Current Source = %.8g\nEnd\n\n'
        'Body Force 2\n  Name = "inject_sink"\n  Current Source = %.8g\nEnd' % (j3, j5))
    return SIF_TEMPLATE.format(bodies="\n\n".join(bodies), materials="\n\n".join(materials),
                                body_forces=body_forces, pin_node=pin_node_1idx)


def run(cmd, cwd):
    print("+ %s" % " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=cwd, env=ELMER_ENV, capture_output=True, text=True)
    print((r.stdout or "")[-2000:])
    if r.returncode != 0:
        raise SystemExit("%s failed:\n%s" % (cmd[0], (r.stderr or "")[-2000:]))
    return r


# ---- minimal raw-appended-VTU reader (no meshio/vtk in system python) ----

def read_vtu(path):
    with open(path, "rb") as f:
        raw = f.read()
    head_end = raw.index(b"<AppendedData")
    head = raw[:head_end].decode("utf-8", "replace")
    npoints = int(re.search(r'NumberOfPoints="(\d+)"', head).group(1))
    ncells = int(re.search(r'NumberOfCells="(\d+)"', head).group(1))
    fields = {}
    for tag in re.finditer(r'<DataArray\b([^>]*)/>', head):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', tag.group(1)))
        if "offset" not in attrs:
            continue
        name = attrs.get("Name") or "__points__"
        fields[name] = (attrs["type"], int(attrs.get("NumberOfComponents", 1)), int(attrs["offset"]))
    blob0 = raw.index(b"_", head_end) + 1
    dtype = {"Float64": np.float64, "Float32": np.float32, "Int32": np.int32, "UInt8": np.uint8}

    def arr(name, n):
        typ, ncomp, off = fields[name]
        pos = blob0 + off
        nbytes = int(np.frombuffer(raw, dtype=np.uint32, count=1, offset=pos)[0])
        a = np.frombuffer(raw, dtype=dtype[typ], count=nbytes // np.dtype(dtype[typ]).itemsize,
                           offset=pos + 4)
        a = a.reshape(-1, ncomp) if ncomp > 1 else a
        assert a.shape[0] == n, (name, a.shape, n)
        return a

    return dict(
        V=arr("potential", npoints),
        points=arr("__points__", npoints),
        geomids=arr("GeometryIds", ncells),
        tets=arr("connectivity", ncells * 4).reshape(ncells, 4),
    )


def elem_field(points, tets, V):
    p = points[tets]  # already metres
    J = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=2)
    Jinv = np.linalg.inv(J)
    g = np.empty((len(tets), 4, 3))
    g[:, 1:, :] = Jinv
    g[:, 0, :] = -g[:, 1:, :].sum(axis=1)
    E = -np.einsum("eki,ek->ei", g, V[tets])
    return np.linalg.norm(E, axis=1)


def main():
    if not os.path.exists(os.path.join(ELMER_BIN, "ElmerSolver")):
        raise SystemExit("Elmer not found at %s -- see fem/README.md 'Running Elmer'" % ELMER_BIN)
    if not os.path.exists(MESH_MSH):
        raise SystemExit("%s missing -- run export_results.py first" % MESH_MSH)
    d = np.load(os.path.join(C.OUT, "solution.npz"), allow_pickle=True)
    nodes, tets, lab = d["nodes"], d["tets"], d["label"]
    order = [str(x) for x in d["order"]]

    def vol_of(name):
        p = nodes[tets[lab == order.index(name)]] * 1e-3
        Jm = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=2)
        return np.abs(np.linalg.det(Jm) / 6.0).sum()

    vol3 = vol_of("contact%d" % C.SOURCE_CONTACT)
    vol5 = vol_of("contact%d" % C.SINK_CONTACT)
    pin_1idx = int(d["pinned_node"]) + 1  # mesh_tagged.msh node ids are 1-indexed

    os.makedirs(WD, exist_ok=True)
    if not os.path.exists(os.path.join(WD, "mesh_tagged.msh")):
        import shutil
        shutil.copy(MESH_MSH, os.path.join(WD, "mesh_tagged.msh"))
    if not os.path.isdir(MESH_DIR):
        run(["ElmerGrid", "14", "2", "mesh_tagged.msh"], cwd=WD)

    sif_path = os.path.join(WD, "case.sif")
    with open(sif_path, "w") as f:
        f.write(build_sif(vol3, vol5, pin_1idx))
    with open(os.path.join(WD, "ELMERSOLVER_STARTINFO"), "w") as f:
        f.write("case.sif\n")
    os.makedirs(os.path.join(WD, "results"), exist_ok=True)
    run(["ElmerSolver", "case.sif"], cwd=WD)

    vtu = os.path.join(WD, "results", "case_t0001.vtu")
    ev = read_vtu(vtu)
    Ve, points, geomids, etets = ev["V"], ev["points"], ev["geomids"], ev["tets"]
    Em = elem_field(points, etets, Ve)
    elab = geomids - 1  # 1-indexed body tag -> 0-indexed ORDER index

    Vp = d["V"].astype(np.float64)
    rng = Vp.max() - Vp.min()
    diff = np.abs(Ve - Vp)
    print("\nCROSS-CHECK  P1 (this repo) vs Elmer StatCurrentSolver, same mesh, same gauge node")
    print("  nodes compared      %d" % len(Vp))
    print("  max |dV|            %.4e V  (%.4f %% of the %.3f V range)"
          % (diff.max(), 100 * diff.max() / rng, rng))
    print("  rms |dV|            %.4e V  (%.5f %%)"
          % (np.sqrt((diff ** 2).mean()), 100 * np.sqrt((diff ** 2).mean()) / rng))
    print("  VERDICT:", "AGREE" if diff.max() / rng < 0.01 else "DISAGREE")

    print("\nElmer field against Khadka 2020 (bipolar 1 A):")
    for name, paper in (("white", 12.0), ("grey", 4.2)):
        i = order.index(name)
        mine = np.percentile(Em[elab == i] / 1e3, 99.9)
        print("  peak |E| %-6s  paper %5.1f kV/m   here %6.2f kV/m (p99.9)  ratio %.2f"
              % (name, paper, mine, mine / paper))
    print("  peak surface voltage  paper 1.2 kV   here %.2f kV (max |V|)"
          % (max(abs(Ve.max()), abs(Ve.min())) / 1e3))


if __name__ == "__main__":
    main()
