# FreeCAD-side FEM route — spinal cord stimulation

The visualizer and `run_all.sh` now use **Elmer StatCurrentSolver with
MPI, Hypre/BoomerAMG and conjugate gradients**. The dedicated **SCS Job Manager**
configures future local jobs: mesh maximum threads, physical-core-capped Elmer
MPI processes and solve memory budget. Older open panels may still show legacy
local controls. These do not read FreeCAD's native Elmer preferences.
With the HPC profile selected, mesh and solve jobs run on the cluster through
SLURM and their results are fetched back.
See [Job Manager](../docs/job_manager.md) and [TODO.md](TODO.md).

Element preparation, sparse reference assembly and independent equation
checks use bounded batches. Production solves avoid building a second global
Python stiffness matrix. The NumPy/SciPy/PyAMG implementation remains a test
reference, not an alternative backend in the visualizer or batch workflow.

See [ELMER_BACKEND.md](ELMER_BACKEND.md) for setup, resource controls,
validation and the preserved pre-Hypre installation.

## Current poster status — September 30, 2026

Both two-lead/16-contact poster models completed local background-inclusive
Elmer solves: dorsal **17,334,566 tets / 9,497.2 s**, ventral
**17,322,431 tets / 9,073.4 s**, each with 14 MPI processes. Mesh workflows took
991.0 and 902.5 s respectively. See the [poster workflow](../docs/poster_fast_workflow.md)
for run IDs, actual mesh controls, counts, residuals and matched-scale figures.
The GUI shows potential in volts; E-field requires gradient post-processing.
Static matched figures are the immediate priority; animation and NEURON results
are not required for this poster.

The roughly 5.6 M classified tets are anatomy/lead only. The other roughly
11.7 M background tets remain in the volume mesh and are included in these
solves. Background is uniform 0.04 S/m, not segmented surrounding tissue.
Unexpected internal unassigned regions have not yet been ruled out. Spatial
classification verification, mesh convergence and background/domain sensitivity
are deferred until after the poster, explicitly recorded in [TODO.md](TODO.md).
Historical validation below must not be generalized to these full models.

## Current workflow and reliability fixes

### Saved FreeCAD documents

Install the lightweight `fem/freecad_startup` directory as a symlink named
`SCSModeling` in FreeCAD's user `Mod` directory (for this installation,
`~/.local/share/FreeCAD/v26-3/Mod/`). Its `Init.py` makes the montage module
loaded before document restoration; `InitGui.py` registers repair callbacks
for older files. Set `SCS_MODELING_REPO` if the checkout moves. The FEM virtual
environment must remain installed with the same Python major/minor as FreeCAD.
Path discovery lives in the imported `scs_bootstrap` module: FreeCAD's addon
loader can execute `Init.py` with its own `__file__`. Loading the callback module
at addon startup also satisfies FreeCAD's document-restore import policy without
changing the trusted-directory list. `test_freecad_startup.py` covers this loader
behavior. A GUI save-copy/reopen check verifies callback and cache restoration.

The generated per-tissue preview surfaces can contain non-manifold edges and
produce “The mesh data structure has some defects” when FreeCAD restores them.
These warnings concern the display surfaces, not a check of the external volume
mesh or Elmer solution. Do not repair these surfaces automatically or interpret
the absence of a solve error as a volume-mesh quality assessment.

Normal File → Save persists currents, the solution path, and mesh cache path/
parameters on the mesh-preview group. Closing and reopening the document then
allows tree current edits followed by Recompute without opening a macro or
running Elmer. The macros target the active document, including suffixed models.
Mesh settings come from that document rather than another model's last run.

Numerical arrays remain in `fem/out/lead_runs/<run>/`; they are **not embedded
in the FCStd file**. Keep that directory alongside the checkout. A missing
solution produces an explicit missing-cache message. This is persistence on
the configured machine, not a self-contained portable FEM archive.

Reuse still checks mesh/solution manifests and current labelled lead geometry.
The comparison tolerates vertex renumbering, binary-STL rounding up to 0.00002 mm,
and alternate triangulations of planar patches with identical boundaries.
It rejects changed mesh settings, missing faces, nonplanar retriangulation, and
meaningful lead movement. Original numerical cache manifests are not rewritten.

Validation: `test_document_cache.py` and `test_montage_persistence.py`, plus a
GUI save-copy/open/recompute test. The reopened object restored both proxies,
matched the existing mesh, and reversed potential signs exactly after reversing
currents, without modifying the solution NPZ.

### Numerical workflow

**The numerical tables below are historical benchmarks.** The default mesher
now covers the union of all anatomical STLs and the selected lead, plus the
configured margin. It no longer reproduces the old canal-only field. This
larger model still needs a mesh-convergence study and independent solver checks
before its fields are used for axon thresholds.

- **SCS Mesh Generator:** mesh sizes, resource estimate, Generate/Abort Mesh,
  tissue classification, dura-interface check and per-tissue preview.
- **SCS Potential Visualizer:** requires a completed matching mesh; owns
  Solve / Plot Field, background-tissue selection, cancellation and montage
  currents. It does not launch a mesher. Changing only currents reuses the
  accepted per-contact basis fields.
- The local solve CPU allocation is capped at available **physical cores** (14
  on the i9-13900H), respecting process affinity. Saved higher values are clamped.
- Closing/reopening the visualizer preserves `SCS_Montage` currents, colour
  clamp, background selection and verified solution path. With the panel closed,
  edit its `ContactN_mA` properties in the tree, then recompute the object to
  update the plots. No Elmer solve is launched by recompute. Geometry or solution
  provenance mismatches are rejected; changing the background scenario can
  require a separate solve. The cached NPZ and run inputs must remain on disk.
- Potential plots include the selected run's contact and insulator surfaces,
  alongside all eleven anatomical classes. Contacts and insulation are grouped
  under Lead 1, Lead 2, … in the potential-field tree. Those lead materials
  were already in the FEM solve. Adding these surfaces changes only plotting;
  it does not invalidate a completed solution. The bulk-tissue colour clamp
  also applies to the lead; `V_volts` retains the unclamped values.
- Choose the potential legend style manually in FreeCAD's legend options.
  The macro does not set the style automatically.
- Size regions include all dura/sheaths, white/grey matter, epidural/CSF,
  roots/DRG/vessels/sympathetic structures, bone/discs and lead hardware.
  `H_Min` and `H_Max` clamp every region's target. They are in **mm**;
  Ansys's displayed metres and its patch-independent refinement settings are
  not interchangeable with these distance-field controls.
- Run keys include anatomy, tissue map, mesh/classification code, dependency requirements,
  lead geometry and mesh parameters. Separate SHA-256 manifests validate
  completed artifacts. Old unmarked runs remain on disk but are not cache hits.
  Lead STLs are stored in each run's own `lead/` directory.
  Solver-code changes invalidate solutions without forcing a remesh.
- Each worker preserves its stdout/stderr in `<run_dir>/<script>.log`, including
  exit status. Silence is logged without killing a potentially healthy worker.
  Abort remains available. An interrupted preview can resume from a matching
  completed volume mesh.

### Memory estimates and interactive display

Meshing uses 80% of currently available RAM as its safety budget; the Job
Manager's explicit memory setting controls the solve, not that mesh gate.
The mesh estimate reports both available RAM and the usable budget and refreshes
while the panel is visible. Closed documents release unused in-memory basis
caches; saved NPZ solutions are retained. A memory estimate is not an OOM guarantee.

Dense mesh and field picking can slow rotation/hover. Disable viewport picking
for these display objects when inspecting dense results; tree visibility still
works. Plot visibility does not change the numerical domain.

### Historical: the 85% abort

The September 15 10:38 kernel log records an OOM kill of the Python worker
at 27,819,292 KiB anonymous resident memory. The saved mesh has **7,630,740
nodes and 47,511,833 tetrahedra**. At 85%, the worker had finished classification
and was about to inspect dura interfaces. The previous code materialized all
tet faces, several sorted copies, and coordinates of all interior faces at once.
Face matching is now partitioned on disk with bounded sort memory, and only
relevant material interfaces need coordinate arrays. Classification also forms
centroids in chunks. Temporary face partitions use `fem/out/`, because `/tmp`
is RAM-backed on this workstation.

Resource estimates now cover mesh **and preview**, compare against currently
available RAM, and are checked again in the worker before meshing. A separate
solve-memory check prevents starting the substantially larger FEM assembly on
an oversized mesh. Full-anatomy count/RAM estimates are provisional, not a
guarantee that every parameter set fits.

### The 60% crash during edge meshing

The September 15 18:38 failure was a segmentation fault during Gmsh 4.15.2
edge meshing, without a corresponding OOM kill. It was reproduced with the
saved Structured size field and 20 boundary-meshing threads (heap corruption).
The mesher now uses one thread for edges and surfaces, retaining the selected
thread count for HXT volume meshing. With those settings, a large reproduction
using that same field completed with 18.65 million tetrahedra and 20 HXT threads.
The exact native race has not been isolated upstream; this is a tested workaround.

`fem/tests/test_gmsh_structured.py` exercises a real Structured field and HXT in
an isolated process and checks positive element volumes and total box volume:

```bash
fem/.venv/bin/python -m unittest discover -s fem/tests -p test_gmsh_structured.py -v
```

### Result acceptance and scenario selection

The solver checks the **true** relative residual for every contact and rejects
nonfinite fields or residuals above **1e-6**. This is an explicit acceptance
tolerance, separate from the requested CG target of **1e-11**. Both tolerances,
contact IDs and measured residuals are stored in the solution; accepted fields
are stored as float64. Disconnected active meshes fail with an actionable error
instead of using one gauge for several independent components. The analytic
verification command fails if **any** resolution fails.
AMG setup uses a fixed random seed. CG continues in batches of 400 iterations
when necessary, up to 1600, without relaxing the acceptance tolerance.

`SCS_BACKGROUND=1 bash fem/run_all.sh` now consistently reads/writes
`solution_bg.npz`, `solution_bg.vtu`, `mesh_tagged_bg.msh`, `field_grid_bg.npz`
and `voltage_slices_bg.png`. It does not export the default scenario by mistake.
Structured-grid and VTU exports support the actual number of basis fields.

Run focused checks with:

```bash
fem/.venv/bin/python -m unittest discover -s fem/tests -v
fem/.venv/bin/python fem/scripts/verify_solver.py
```

Still outstanding: anisotropic white matter, the ambiguous root “Middle” tissue
assignment, anatomical resolution/convergence checks and NEURON coupling.
No conductivity was changed in this reliability work.

## Historical canal-only workflow and validation

An independent second opinion on the RADO-SCS current-flow problem, built with
open tools only (gmsh + a P1 finite-element solver in this repo), deliberately
sharing nothing with the Ansys pipeline except the geometry and the
conductivities in `src/ansys/tissue_map.yaml`. If this and Ansys converge on a
similar field, that is mutual validation; if they do not, the disagreement is
informative.

Everything in this directory is generated by scripts. Nothing was done by GUI
clicking, and nothing outside `fem/` was modified.

## Setup

`requirements.txt`'s packages (scipy, gmsh, pyamg, meshio) have no wheels yet
for the system Python (3.14 as of 2026-09), and `pip install --user` is
refused outright (PEP 668, "externally managed environment") since this is
Ubuntu's own Python, not a user one. Needs its own venv, once:

    python3 -m venv fem/.venv
    fem/.venv/bin/pip install -r fem/requirements.txt

`fem/run_all.sh` picks up `fem/.venv` automatically if it exists (set `PY=`
to override). No `sudo` needed anywhere in this.

    bash fem/run_all.sh

**Result: a voltage field, solved and cross-checked. Peak |E| in white matter
came out at 11.5 kV/m against Khadka's published 12 kV/m on the five-compartment
model; 10.0 kV/m since all 240 anatomical bodies were added on 2026-09-15. See
RESULTS, which opens with what that change moved.**

![voltage on a sagittal and axial slice through the lead](out/voltage_slices.png)

## What the equation is

    ∇·(σ∇V) = 0 ,      J = −σ∇V

Steady current conduction. Solved directly, in volts, with linear tetrahedra.
No physics analogy is used for the field of record.

## Solver choice, and why

**Chosen: a P1 (linear-tetrahedron) FEM written in `fem/scripts/`, solved with
scipy + pyamg.**

The options and how they scored:

| option                                             | verdict                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| -------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Elmer `StatCurrentSolver`**                       | The _right_ answer, run directly (not through FreeCAD's GUI) via `fem/scripts/crosscheck_elmer.py`: `ElmerGrid` converts `mesh_tagged.msh` to Elmer's native format, then `ElmerSolver` solves the real electric-conduction equation with no analogy at all — `Electric Conductivity` per body in S/m, a volumetric `Current Source` in A/m³ on the driven contacts, the same single-node gauge pin as the P1 solve. **Run 2026-09-11: agrees with the P1 field to 9.1e-3 V max, 0.0003 % of the 2667.6 V range.** See *Running Elmer* below and RESULTS.                                                     |
| **CalculiX steady-state heat as an exact analogy** | Sound physics — ∇·(k∇T)=0 _is_ ∇·(σ∇V)=0 — and `ccx` is installed. Used here, but as a **cross-check only** (`fem/scripts/crosscheck_ccx.py`), never as the field of record. The reason is exactly the trap named in the brief: a `.frd` whose nodal field is called `NDTEMP` but means volts is a landmine for whoever reads it next. Note that it is _not_ rejected for the floating contacts — the high-conductivity-body treatment used here needs no special constraint and works in any solver, which the cross-check demonstrates by reproducing the same field. |
| **P1 FEM in this repo**                            | Chosen. It solves the equation literally, in volts; it makes the floating contacts exact rather than approximate; it can hand the NEURON stage V at arbitrary points without a format round-trip; and — the actual point of this exercise — it shares no code with Ansys, so agreement between the two means something. Its risk is that it is new code, so it is verified against a closed-form solution before it is believed (see _Verification_).                                                                                                                   |

### Running Elmer

Installed at **`/opt/elmerfem`**, version **26.2**, built from source (the
`elmer-csc` PPA does not cover this Ubuntu release). It is NOT on `PATH` and
nothing has been added to a shell profile, so every session needs:

    export PATH=/opt/elmerfem/bin:$PATH
    export ELMER_HOME=/opt/elmerfem
    export ELMER_SOLVER_HOME=/opt/elmerfem/share/elmersolver

**Correction to what was believed here before actually running it: `PATH`
alone is NOT sufficient.** Without `ELMER_SOLVER_HOME`, `ElmerSolver` fails
immediately with `InitializeElementDescriptions: elements.def not found` —
it does not fall back to a path relative to the binary. All three variables
above are needed; `fem/scripts/crosscheck_elmer.py` sets them itself so no
shell setup is required to reproduce the run. Present and relevant:

    /opt/elmerfem/bin/ElmerSolver          serial (symlinked to ElmerSolver_mpi)
    /opt/elmerfem/bin/ElmerGrid            mesh conversion
    /opt/elmerfem/share/elmersolver/lib/StatCurrentSolve.so      <- the solver used
    /opt/elmerfem/share/elmersolver/lib/elements.def             <- needs ELMER_SOLVER_HOME to be found

FreeCAD's FEM workbench will not find Elmer by itself: point
**Edit → Preferences → FEM → Elmer** at `/opt/elmerfem/bin/ElmerSolver` if
driving it from the GUI. This repo does not — see below.

`fem/out/mesh_tagged.msh` is written with one physical volume per tissue
specifically so it can go straight into `ElmerGrid 14 2 mesh_tagged.msh`,
which is exactly what `crosscheck_elmer.py` does before writing a `.sif` and
running `ElmerSolver` — the **third** independent solve on the identical
mesh. `ElmerGrid` reported "No lower dimensional elements present!" because
`mesh_tagged.msh` carries only 3D physical volumes, no 2D physical surfaces —
so `mesh.boundary` comes out empty and there is nothing to attach a surface
BC to. That is not a problem here: the outer canal wall is insulating by
Elmer's default (zero-flux Neumann, the same "do nothing" boundary the P1
solve relies on), the driven contacts get a volumetric `Current Source`
instead of a surface condition, and the gauge pin targets a node directly
(`Target Nodes`) rather than a boundary. MUMPS (linked into this build) is
used as a direct solver rather than iterative CG/ILU, given the ~5×10⁸
conductivity contrast between the clamped metal and the lead insulation —
Elmer's own timer reports 79 s for the solve (CPU time summed across MUMPS's
internal threads), but wall time for the whole script — mesh conversion,
assembly and solve — was 20 s. Reproduce with:

    python fem/scripts/crosscheck_elmer.py

## Geometry, and what fought

### The good news

All six bodies of the stripped model are clean. Measured with
`fem/scripts/stlio.py`:

| body                          | triangles    | closed | non-manifold edges | shells                             | volume (mm³)         |
| ----------------------------- | ------------ | ------ | ------------------ | ---------------------------------- | -------------------- |
| `neuro_EpiduralSpace-1`       | 3512         | yes    | 0                  | 1                                  | 11295.4              |
| `neuro_Meninges-1` (dura)     | 3560         | yes    | 0                  | 1                                  | 2151.0               |
| `neuro_CSF-1`                 | 4276         | yes    | 0                  | 1                                  | 8835.8               |
| `neuro_whitemater-1`          | 10770        | yes    | 0                  | 1                                  | 2295.3               |
| `neuro_GreyMater-1`           | 8354         | yes    | 0                  | 1                                  | 995.6                |
| lead (8 contacts + insulator) | 8×500 + 4500 | yes    | 0                  | 1 per contact, 9 for the insulator | 3.68–3.71 each, 23.4 |

No non-manifold edges (every edge shared by exactly two facets), no boundary
edges, no zero-area facets, and a divergence-theorem volume that agrees with a
Monte-Carlo occupancy integral to under 1 % — which a self-intersecting or
inconsistently wound shell would not manage. Self-intersection was not tested
directly. The STL quality worry in the brief did not materialise for this
subset. Every body in `NBF_RADO-SCS_fem.FCStd` has an
identity `Placement`, so the STL files in `STL_files/` are already in document
coordinates and the generated lead drops in untransformed, as promised.

### The actual blocker

**The shared interfaces are not tessellated the same way.** Counting welded
vertices that coincide to 1e-4 mm between neighbouring bodies:

    epidural ∩ dura      4  shared vertices (of ~1760)
    dura     ∩ csf       4
    csf      ∩ white   237  (of 2138)
    white    ∩ grey    442

The compartments tile space geometrically — confirmed by Monte-Carlo
integration: sampled occupancies reproduce the divergence-theorem volumes to
better than 1 %, and pairwise overlap totals 20 mm³ out of 25 573 mm³, i.e.
0.08 % (epidural∩dura 16, dura∩CSF 2.3, CSF∩white 1.5, white∩grey 0.5 mm³),
all of it tessellation slivers at interfaces — **but their surface meshes do not
match node-for-node.** There is therefore no way to assemble them into a
conformal multi-volume mesh without a boolean/remesh repair step. This is
almost certainly the same wall the Ansys side keeps hitting.

### What was done instead

One conformal tetrahedral mesh of a single box containing the whole canal,
graded by distance to the structures that matter, with **a tissue assigned to
every tetrahedron by point-in-shell classification** of its centroid
(`fem/scripts/inside.py`, ray casting with an XY grid index; validated against
the divergence-theorem volumes above). Tets that land outside every body are
deleted, which leaves the tissue union with a naturally insulating outer
surface.

Consequence, stated plainly: **interfaces are resolved to the local element size
rather than followed exactly.** That is the headline approximation of this
route. It is why the mesh is refined hardest on the dura.

### The dura is the expensive part

Measured along a dorsoventral ray through the centre of contact 3:

    epidural fat   2.35 mm     ← the lead sits in here
    dura           0.50 mm
    CSF            2.65 mm
    white matter   2.45 mm
    grey matter    0.35 mm  (ray clips the edge of the horn at this x)

The dura is 0.5 mm thick and, at σ = 0.037 S/m, is the most resistive tissue in
the model and the barrier the current has to cross. Its mean thickness over the
whole body (2V/A) is 0.48 mm, so 0.5 mm is representative, not a local thin
spot. Resolving it is what sets the mesh size, and a mesh too coarse there
does not merely blur the answer — it can punch holes through the barrier and
short epidural fat straight to CSF. `assign_and_solve.py` therefore measures
that leakage explicitly rather than assuming it away (see RESULTS).

## Boundary conditions, as actually imposed

- **+1 A into contact 3, −1 A out of contact 5**, numbered as the STL filenames
  run. That numbering runs caudal → rostral: contact 1 spans z = 123.25–126.41 mm
  and contact 8 spans z = 151.32–154.34 mm, so **contact 3 (z = 131.27–134.39) is
  the caudal member of the driven pair and contact 5 (z = 139.29–142.37) the
  rostral one**, 8.0 mm apart centre to centre, with contact 4 floating between
  them. Injected as a volumetric source spread over the nodes of the contact
  body, weighted by nodal volume. Because a contact is metal, it is an
  equipotential body and the distribution inside it does not affect the exterior
  field; the metal puts the current on its own surface. The equipotentiality is
  measured, not assumed.
- **The six undriven contacts are exact floating conductors.** They are left as
  ordinary high-conductivity bodies with no source term. With no source inside,
  conservation forces the net current through each to be zero — that is exact,
  not an approximation — and the high conductivity makes each equipotential.
  Nothing was silently dropped. Both properties are reported per contact.
- **Conductivity of the metal is clamped to 1×10⁴ S/m**, not the 4×10⁶ S/m in
  `tissue_map.yaml`. A 10⁸ conductivity ratio against epidural fat wrecks the
  conditioning of the linear system for no physical gain: at 10⁴ S/m a contact
  is already equipotential to ~10⁻⁵ of the driving voltage, which the reported
  spread confirms. `config.SIGMA_METAL_TRUE` keeps the real value so the clamp
  is auditable. Nothing else departs from `tissue_map.yaml`.
- **Outer boundary insulating.** Zero normal current is the natural ("do
  nothing") condition of the weak form, so it needs no code: every face with no
  neighbouring element carries zero current. The boundary sits at the outer
  surface of the epidural body, i.e. at the canal wall.
- **The gauge.** A pure all-Neumann problem is singular up to an additive
  constant. One node — the mesh node furthest from the lead axis — is pinned to
  0 V. Because the sources sum to exactly zero, that node draws no real
  current; the reaction there is reported as a fraction of 1 A, and is the check
  that the pin is a gauge choice rather than a current path.

### Basis fields, for the NEURON stage

Eight solves are done, one per contact: +1 A into contact _i_, −1 A spread over
the outer boundary weighted by nodal area. Any zero-net-current montage is then
a linear combination, V = Σ Iᵢ φᵢ, in which the boundary return cancels
identically. The requested bipolar case is φ₃ − φ₅ at 1 A. This matches the
shape of the contract in `docs/neuron_plan.md` chunk 3 (one φ per contact,
explicit origin / spacing / shape / units) without pre-empting its unit
convention: everything here is mm and volts, and the NEURON side converts once,
to µm and mV, in its own `field.py`.

## Verification

`fem/scripts/verify_solver.py` solves a spherical shell a ≤ r ≤ b with current I
injected at r = a and r = b held at 0, against the closed form

    V(r) = I / (4πσ) · (1/r − 1/b)

at three mesh densities. This exercises assembly, the mm → m unit conversion,
the volumetric injection and the Dirichlet pin together.

    h at r=a   nodes    L∞ error   RMS error   V(a) FEM / exact
    0.60 mm     5423     2.07 %     0.66 %     54.24 / 55.26 V
    0.40 mm    16048     1.14 %     0.33 %     54.77 / 55.26 V
    0.28 mm    41599     0.61 %     0.17 %     55.01 / 55.26 V

Converging, with the correct absolute voltage — so the units are right, not just
the shape. Current in = +1.000000 A, out = −1.000000 A at every density.

**This test earned its keep.** The first version of the solver transposed the
inverse Jacobian when forming the basis-function gradients. The resulting matrix
was still symmetric with zero row sums, so **the current balance still came out
exactly ±1.000000 A while the field was wrong by a factor of 2.2**, and the
error did not shrink under refinement. A conservation check alone would have
passed a wrong answer. The comment at that line in `assign_and_solve.py` records
this.

## RESULTS

A voltage field was obtained. Everything below is measured output from
`fem/out/`, not an estimate.

### 2026-09-15 — the model is no longer five compartments, and the numbers moved

Everything in the rest of RESULTS describes the **five-compartment** model
(epidural, dura, CSF, white, grey) that `fem/out/solution.npz` still holds, and
that the CalculiX and Elmer cross-checks were run against. It is kept because
that is what those committed files are.

`config.TISSUE_BODIES` now classifies **all 240 anatomical STLs** in
`STL_files/` into the eleven tissue classes `src/ansys/tissue_map.yaml` defines
— vertebrae, discs, vasculature, nerve roots, DRG and the sympathetic chain
included — each with its own conductivity read from that file. See
*What is approximate* item 3. `bash fem/run_all.sh` therefore no longer
reproduces the table below; it reproduces this one. Both were produced by the
identical `assign_and_solve.main()`, on the identical `fem/out/mesh.npz`, and
differ only in how many bodies the tets were classified against:

| quantity                        | 5 compartments | all 11 tissue classes | change  |
| ------------------------------- | -------------- | --------------------- | ------- |
| active nodes                    | 289 136        | 301 693               | +4.3 %  |
| active tets                     | 1 738 750      | 1 825 346             | +5.0 %  |
| tets dropped as background      | 148 943        | 62 347                | −58 %   |
| dura leak                       | 0.000 %        | **0.000 %**           | —       |
| bipolar transfer impedance      | 2667.6 Ω       | **1737.1 Ω**          | −34.9 % |
| peak \|E\| white matter (p99.9) | 11.51 kV/m     | 10.02 kV/m            | −12.9 % |
| peak \|E\| grey matter (p99.9)  | 7.91 kV/m      | 7.04 kV/m             | −11.0 % |
| peak surface voltage            | 1.37 kV        | 0.90 kV               | −34 %   |
| classification noise            | 0.0021 %       | 0.0608 %              | ×29     |
| field-grid in-mesh coverage     | 52.4 %         | 79.9 %                | —       |

Tets by tissue, and what that says about whether a structure is really there:

| tissue            | tets    | mm³      | note                                   |
| ----------------- | ------- | -------- | -------------------------------------- |
| epidural          | 468 782 | 8 958.3  |                                        |
| dura              | 622 299 | 5 543.1  | main meninges **+ 24 root sheaths**     |
| CSF               | 512 910 | 8 384.1  | main CSF + 8 DRG coatings              |
| white             | 69 258  | 2 294.0  | unchanged                              |
| grey              | 30 154  | 997.9    | unchanged                              |
| vertebra          | 60 311  | 34 295.5 | was background; drives most of the −35 % |
| root              | 22 651  | 496.8    |                                        |
| blood             | 18 481  | 333.8    | radicular arteries, 0.66 S/m            |
| disc              | 12 884  | 8 938.8  |                                        |
| DRG               | **126** | 11.3     | **barely resolved** — see below         |
| sympathetic chain | **0**   | 0        | **entirely outside the mesh box**       |

The impedance drop is the expected direction and is larger than the ±16 % the
`SCS_BACKGROUND=1` variant measured, for the same reason and more of it: the
return current now has 34 295 mm³ of vertebral bone at 0.04 S/m, plus roots,
discs and vessels, to spread into instead of stopping at the canal wall. The
two published quantities move apart from Khadka rather than towards him (white
matter 0.96 → 0.84 of the paper, surface voltage 1.14 → 0.75); that is not
evidence the richer model is worse, because this is still isotropic white
matter and a very coarsely resolved periphery — but it is the honest reading
and it should not be dressed up. The CalculiX and Elmer cross-checks have
**not** been re-run against this model.

**The DRG and the sympathetic chain are present in name only.** The mesh box is
the epidural envelope + 1 mm, so the sympathetic chain lies entirely outside it
and gets zero tets; and the size field still grades on only four regions (dura,
lead, white, canal), so a DRG that occupies 11 mm³ collects 126 tets and has
5.6 % of them isolated from every face neighbour — noise, not anatomy. See
*What is approximate* items 3 and 10.

### The mesh

gmsh 4.15 (HXT parallel Delaunay), one box graded by four distance fields,
**308 081 nodes / 1 887 693 tetrahedra in 760 s** — of which 12 min is the
distance-field precomputation (cached in `fem/out/sizefield.npz`, so a re-run is
~2 min) and 2.1 s is the actual 3D meshing. Deleting the 148 943 tets that fall
outside every tissue leaves the solved domain:

    289 136 nodes,  1 738 750 tetrahedra

Target element size: 0.25 mm on the dura, 0.30 mm on the lead, 0.50 mm in the
cord, 0.70 mm at the canal wall, 2.5 mm in the discarded surround.

### Did the classification actually reproduce the anatomy?

Tet-sum volume against the STL's own divergence-theorem volume:

| tissue     | tets         | meshed vol mm³ | STL vol mm³    | error   |
| ---------- | ------------ | -------------- | -------------- | ------- |
| epidural   | 584 563      | 11 247.34      | 11 295.4       | −0.43 % |
| dura       | 519 027      | 2 150.63       | 2 151.0        | −0.02 % |
| CSF        | 528 258      | 8 829.16       | 8 835.8        | −0.08 % |
| white      | 69 258       | 2 294.02       | 2 295.3        | −0.06 % |
| grey       | 30 154       | 997.94         | 995.6          | +0.23 % |
| 8 contacts | 498–545 each | 3.67–3.79 each | 3.68–3.71 each | ≤ +2 %  |
| insulator  | 3 335        | 23.55          | 23.4           | +0.6 %  |

**The dura barrier leaks nowhere.** Of 13 893.7 mm² of resolved dura interface,
the area of faces putting epidural fat or the lead directly against CSF, white
or grey matter — a short across the most resistive tissue in the model — is
**0.000 mm², 0.000 %**. This was the single biggest risk in the approach and it
did not materialise; 519 k tets in a 0.5 mm shell is enough.

**Classification noise is negligible.** `inside.py` decides with one +Z ray, and
the lead-designer agent documented that single-ray parity can invert on a
grazing edge. Counting tetrahedra whose tissue differs from _every_ face
neighbour — the signature of a miscount rather than anatomy — gives **37 out of
1 738 750, or 0.0021 %**, the worst single tissue being white matter at
0.023 %. Isolated elements at that rate cannot move a diffusion solution.

### The solve

Eight basis solves (one per contact), CG preconditioned by smoothed-aggregation
AMG, ~35 s each, 300 s total including assembly.

- Relative residual **2.6×10⁻⁸ to 1.4×10⁻⁷**. Every solve hit the 400-iteration
  cap rather than the 10⁻¹¹ tolerance — the 5×10⁸ spread in σ (2×10⁻⁵ for the
  lead insulation up to 10⁴ for metal) makes the system stiff. 10⁻⁷ is far
  tighter than the discretisation error, so this is not a limitation on the
  answer, but it is not full convergence and is stated as such.
- **Gauge pin reaction: ≤ 2.4×10⁻¹⁰ A** against a 1 A drive, across all eight
  solves. The pin is a gauge choice, not a current path, as intended.
- Global sum of nodal currents 6.4×10⁻¹² A.

### Boundary conditions, verified rather than asserted

| contact | role       | potential      | equipotential spread | net current     |
| ------- | ---------- | -------------- | -------------------- | --------------- |
| 1       | floating   | +34.10 V       | 0.0002 % of drive    | +0.000000 A     |
| 2       | floating   | +204.29 V      | 0.0010 %             | +0.000000 A     |
| **3**   | **source** | **+1297.67 V** | 0.0028 %             | **+1.000000 A** |
| 4       | floating   | −31.79 V       | 0.0032 %             | +0.000000 A     |
| **5**   | **sink**   | **−1369.91 V** | 0.0023 %             | **−1.000000 A** |
| 6       | floating   | −239.48 V      | 0.0010 %             | +0.000000 A     |
| 7       | floating   | −80.84 V       | 0.0002 %             | +0.000000 A     |
| 8       | floating   | −57.16 V       | 0.0000 %             | +0.000000 A     |

The six floating contacts are equipotential to better than 0.004 % of the
driving voltage. The 10⁴ S/m clamp is therefore doing its job and the 4×10⁶ S/m
of the real metal would buy nothing.

Read the net-current column for what it is. Zero net current through a floating
contact is _imposed by construction_ here — those bodies simply have no source
term — so the column is a convergence check on the linear solve, not an
independent discovery. The equipotential-spread column is the one that measures
something that could have come out badly, and it is the reason the clamp can be
called adequate. The genuinely exact part of the floating treatment is the
physics, not the number: with no source inside a body, conservation leaves it
no choice.

**Bipolar transfer impedance: 2667.6 Ω** between contacts 3 and 5.

### The field, against the one published number we can check

Khadka et al. 2020 report, for the full 19-compartment model with anisotropic
white matter at bipolar 1 A:

| quantity                   | Khadka 2020 | this stripped model    | ratio |
| -------------------------- | ----------- | ---------------------- | ----- |
| peak \|E\| in white matter | 12 kV/m     | **11.51 kV/m** (p99.9) | 0.96  |
| peak \|E\| in grey matter  | 4.2 kV/m    | 7.91 kV/m (p99.9)      | 1.88  |
| peak surface voltage       | 1.2 kV      | 1.37 kV                | 1.14  |

White matter and the electrode voltage land within 4 % and 14 % of the paper.
That is closer than this model deserves and should not be read as validation of
the details — it is a stripped 6-compartment model with isotropic white matter
against a 19-compartment one at 150 M elements — but it does say the geometry,
the conductivities, the current injection and the units are all in the right
place. Grey matter is high by 1.9×, which is the direction isotropic white
matter would push it: with no longitudinal shunt at 0.6 S/m, more current
crosses into the grey horns instead of running along the columns.

Percentiles are quoted rather than maxima because the element maximum sits on
the singular edge of a contact: `|E|` maxima reach 3.0×10³ kV/m in the epidural
fat and 2.4×10³ kV/m in the dura, both immediately adjacent to electrode
corners, where a sharp-edged conductor in a P1 mesh has no finite field. The
p99.9 figures are away from those corners.

### Independent second solver

`fem/scripts/crosscheck_ccx.py` re-solves the identical mesh and boundary
conditions in **CalculiX 2.21** through the steady-heat analogy, iterative
Cholesky, 24.5 s:

    nodes compared   289 136
    max |ΔV|         6.19e-03 V   = 0.0002 % of the 2667.6 V range
    rms |ΔV|         1.42e-04 V   = 0.00001 %
    VERDICT          AGREE

Two independently written FEM codes, same mesh, agreeing to two parts in a
million. Combined with the analytic sphere test above, the solver is not the
weak link in this pipeline; the geometry and the stripped compartment list are.

### Historical independent cross-check — Elmer

`fem/scripts/crosscheck_elmer.py` re-solves the identical mesh in **Elmer
26.2's own `StatCurrentSolver`** — unlike the CalculiX check, this is the real
electric-conduction equation in Elmer's own solver, not a physics analogy.
MUMPS direct solve, 20 s wall time including mesh conversion:

    nodes compared   289 136
    max |ΔV|         9.12e-03 V   = 0.0003 % of the 2667.6 V range
    rms |ΔV|         7.94e-03 V   = 0.0003 %
    VERDICT          AGREE

and, computing the field the same way `analyze.py` does (element-constant
`E = -∇V`, p99.9 to avoid the singular contact edges):

| quantity                   | Khadka 2020 | P1 (this repo) | Elmer      |
| --------------------------- | ----------- | --------------- | ---------- |
| peak \|E\| in white matter | 12 kV/m     | 11.51 kV/m      | 11.51 kV/m |
| peak \|E\| in grey matter  | 4.2 kV/m    | 7.91 kV/m       | 7.91 kV/m  |
| peak surface voltage       | 1.2 kV      | 1.37 kV         | 1.37 kV    |

Three independently written codes — this repo's P1 solver, CalculiX via the
heat analogy, and now Elmer's native electric-conduction solver — land on the
same field to 4-5 significant figures. That is not "roughly similar," it is
the same linear system solved three different ways. It closes the loop
`HANDOFF.md` opened: Elmer was installed specifically to give a solve
with no analogy and no shared code, and it reproduces the P1 field, not just
its ballpark. The disagreement with Khadka's grey-matter number is real
physics (isotropic white matter, see above and the anisotropy item in "What
is approximate"), not a solver artifact — all three solvers agree on it.

### How much does stopping at the canal wall cost?

The insulating boundary sits at the outer surface of the epidural body, so no
return current can spread into bone or paraspinal soft tissue. `SCS_BACKGROUND=1`
reruns the identical mesh with the discarded surround kept and filled at the
`tissue_map.yaml` vertebra value of 0.04 S/m (64 839 mm³ of it), which is the
opposite extreme — a fully conductive bone block out to a rectangular box.
The truth is in between. Written to `fem/out/solution_bg.npz`:

| quantity                        | insulating at canal wall | 0.04 S/m surround | change  | Khadka   |
| ------------------------------- | ------------------------ | ----------------- | ------- | -------- |
| bipolar transfer impedance      | 2667.6 Ω                 | 2238.5 Ω          | −16.1 % | —        |
| peak \|E\| white matter (p99.9) | 11.51 kV/m               | 9.65 kV/m         | −16.2 % | 12 kV/m  |
| peak \|E\| grey matter (p99.9)  | 7.91 kV/m                | 6.66 kV/m         | −15.7 % | 4.2 kV/m |
| peak surface voltage            | 1.37 kV                  | 1.12 kV           | −18.2 % | 1.2 kV   |

So the outer boundary treatment is worth about **16 %** on everything, and the
two treatments bracket the paper on both published quantities — surface voltage
(1.12 / 1.37 kV around 1.2) and white-matter field (9.65 / 11.51 around 12). The
field of record is the insulating one, because that is what the brief specified;
the variant is here to bound the choice, not to be averaged with it. Adding real
vertebrae is the way to close this, and it is the obvious next compartment to
un-strip.

### Handover to NEURON

`fem/out/field_grid.npz` — origin, 0.25 mm spacing, shape (93, 77, 243),
`phi[8]` in V/A, tissue label per grid point, units strings, and the montage
note. In-mesh coverage is 52.4 %; the rest of that box is outside the canal and
comes back `NaN`, which is the honest answer. `fem/scripts/sample_example.py`
demonstrates both routes on a straight dorsal-column fibre and cross-checks
them: **grid trilinear vs exact tetrahedral sampling agree to 0.09 V max,
0.023 V rms, over a 61.5 V swing along the fibre (0.15 %)**, and the activating
function peaks at z = 132.0 mm, inside contact 3 (131.27–134.39 mm), which is
where it should be.

**One warning for the NEURON stage, found by running it.** The activating
function is a second difference, and a P1 field is only C⁰ — its second
derivative is a set of jumps at element faces. Near the contacts the signal
(≈3.5×10⁶ V/m²) dominates easily, but 10–20 mm away, where d²V/ds² falls to
≈10⁵ V/m², sampling at 0.5 mm produces visible sign flips between adjacent
points. This is the failure mode `docs/neuron_plan.md` predicts ("if it is
noisy, the FEM grid is too coarse"), and it is not fixed by a finer _sampling_
step — that makes it worse. The fixes, in order of preference: refine the cord
mesh below the current 0.5 mm target; move to quadratic tetrahedra; or fit a
smooth function to V along each trajectory before differencing. Worth settling
in chunk 5, before thresholds are computed on top of it.

## Files

    fem/README.md                  this file
    fem/run_all.sh                 reproduce everything
    fem/requirements.txt           pip deps (no root needed)
    fem/scripts/config.py          paths, the tissue set + classification priority and
                                   every conductivity (all READ from tissue_map.yaml), BC spec
    fem/scripts/stlio.py           binary/ASCII STL reader + topology audit
    fem/scripts/inside.py          point-in-closed-mesh test, XY-grid accelerated
    fem/scripts/build_mesh.py      Gmsh graded tet mesh of the anatomical domain
    fem/scripts/assign_and_solve.py tissue classification and Python regression reference
    fem/scripts/element_batches.py bounded element calculations, reference assembly, matrix action
    fem/scripts/elmer_backend.py   production Elmer/Hypre MPI solver and independent checks
    fem/scripts/solve_resources.py local CPU/memory budget and process-tree monitor
    fem/scripts/solve_lead.py       live-lead worker with provenance validation
    fem/scripts/run_elmer.py        batch-workflow entry for the same Elmer backend
    fem/scripts/field.py           tet-mesh sampling at arbitrary points; grid export
    fem/scripts/export_results.py  VTU + tagged .msh + field_grid.npz
    fem/scripts/analyze.py         per-tissue |E| stats, Khadka comparison, noise
    fem/scripts/sample_example.py  worked handover to the NEURON stage
    fem/scripts/plot_slices.py     the figure
    fem/scripts/verify_solver.py   analytic sphere verification
    fem/scripts/crosscheck_ccx.py    independent CalculiX second opinion (heat analogy)
    fem/scripts/crosscheck_elmer.py  independent Elmer third opinion (real electric-conduction solver)

  interactive (FreeCAD) half -- see each file's own docstring
    fem/scripts/SCS_Mesh_Generator.FCMacro       size-field dock: live cost estimate,
                                                 Generate Mesh and per-tissue preview
    fem/scripts/SCS_Potential_Visualizer.FCMacro montage dock: Plot Field on the live lead
    fem/scripts/live_lead.py       finds/tessellates/fingerprints the lead in the document tree
    fem/scripts/async_runner.py    one QProcess + stall watchdog + progress parser, shared
    fem/scripts/mesh_preview.py    Stage 1 subprocess: build, classify, dura check, preview.vtp
    fem/scripts/solve_lead.py      local Elmer/Hypre worker for the selected lead
    fem/scripts/mesh_cost.py       predicts tets/RAM/time for a parameter set (calibrated)
    fem/scripts/mesh_cost_server.py  keeps one CostModel alive and answers over a pipe
    fem/scripts/scs_montage.py, scs_montage_feature.py  superposition + the document object
    fem/out/voltage_slices.png     THE FIGURE: sagittal + axial through the lead
    fem/out/solution.npz           nodes, tets, tissue, sigma, phi[8], V -- field of record
    fem/out/field_grid.npz         structured phi grid for NEURON (see field.py)
    fem/out/solution.vtu           ParaView / FreeCAD FEM post-processing
    fem/out/mesh_tagged.msh        gmsh 2.2, one physical volume per tissue (for Elmer)
    fem/out/mesh.msh, mesh.npz     the raw mesh
    fem/out/sizefield.npz/.dat     cached size field -- delete to force a rebuild
    fem/out/elmer/                 ElmerGrid mesh + case.sif + results/case_t0001.vtu

## What is approximate, in one list

1. **Interfaces are resolved to element size, not followed exactly** — the
   consequence of classifying tets instead of meshing bodies conformally. The
   dura leakage figure in RESULTS bounds how much this costs where it matters.
2. **White matter is isotropic** at 0.1432 S/m, following `tissue_map.yaml`.
   Khadka's paper uses 0.1432 transverse / 0.6 longitudinal. For dorsal-column
   fibre recruitment this matters and should be added before any comparison with
   the paper's thresholds is taken seriously.
3. **~~Stripped model.~~ FIXED 2026-09-15 — but read the two caveats.** All 240
   anatomical STLs are now classified, into the eleven tissue classes
   `src/ansys/tissue_map.yaml` defines, each with that file's own conductivity
   (`config.TISSUE_BODIES` / `config.ANATOMY_CLASS`). Vertebrae, discs,
   vasculature, nerve roots, DRG and the sympathetic chain are in. What that
   cost and what it changed is in RESULTS, first subsection. Two caveats:
   - **The mesh box is still only the epidural envelope + 1 mm**, so structures
     outside the canal are clipped or absent entirely — the sympathetic chain
     gets **zero** tets. "Classified" is not the same as "meshed". TODO:
     size the box to the union of the bodies actually being classified.
   - **Only four regions drive element SIZE** — see item 10.
4. **TODO, UNRESOLVED, and it sits right next to the electrodes: are the roots'
   "…Middle…" bodies CSF or nerve?** RADO models each root as
   Inside / Middle / OutsidMenging, mirroring the DRG's
   in_middle / csf_coating / menging_coating triple. If that analogy holds, the
   Middle bodies are CSF at 1.7 S/m, not nerve at 0.1432 S/m — a **12×**
   conductivity difference in bodies a few millimetres from the contacts.
   `src/ansys/tissue_map.yaml` raises this question itself and currently assigns
   them `nerve_root`; `fem/scripts/config.py` follows tissue_map.yaml and
   deliberately does **not** invent an answer. Settle it in tissue_map.yaml, by
   measurement or by asking whoever built the CAD, not in the FEM code.
5. **The domain stops at the canal wall** with an insulating boundary. Real
   return current spreads into bone and paraspinal soft tissue. Measured on the
   five-compartment model: this choice is worth about 16 % on impedance and on
   the cord field — see _How much does stopping at the canal wall cost?_ above.
   The two treatments bracket the paper's published numbers, so it is bounded,
   not unknown. Item 3 has now put real vertebrae inside that boundary, which
   moved the impedance by −35 %; the boundary itself has not moved.
6. **Metal conductivity clamped** to 1e4 S/m as described above.
7. **The model is truncated in z** by RADO itself at z = 59.4 and 164.9 mm.
   Contact 8 ends 10.6 mm from the rostral cut and the insulator 4.6 mm from it,
   so the rostral boundary is close enough to the lead to matter for contact 8.
   Contacts 3 and 5, the pair actually driven, are 30+ mm away from either cut.
8. **The activating function is noisy far from the contacts.** A P1 field is
   C⁰, so its second derivative is a set of jumps at element faces. Near the
   contacts this is invisible; 10–20 mm away it is not. See the warning at the
   end of RESULTS, and settle it before thresholds are computed on top of it.
9. **The eight basis solves stopped on an iteration cap**, at 10⁻⁷ relative
   residual rather than the 10⁻¹¹ requested. That is well below the
   discretisation error, but it is not full convergence.
10. **All tissues are classified; only FOUR regions drive element size.**
    `build_mesh.py`'s size field is graded against exactly four hardcoded
    surfaces — dura, lead, white matter, canal wall — and `config.ORDER` /
    `classify_tets()` is a separate mechanism that decides which σ each
    tetrahedron gets. So every tissue added in item 3 is meshed at whatever
    ambient size those four regions happen to produce nearby, with no
    refinement of its own. For a thin or small structure that means
    under-resolution, or being missed outright by centroid classification when
    the local element is bigger than the structure: measured on the
    field-of-record mesh, the DRG gets 126 tets for 11 mm³ with 5.6 % of them
    isolated from every neighbour, the radicular arteries 18 481 tets with
    0.54 % isolated, against 0.02 % for white matter. The radicular arteries
    sit right beside the electrodes at 0.66 S/m, so this is not academic.
    TODO: per-tissue-group size control in tiers, replacing the four fixed
    regions — Mohamed's own Ansys setup already scopes a meshing method per
    Named Selection (Fine/Moderate/Coarse-grained Tissue, SCS Lead, Meninges),
    each with its own max element size, and that is the shape to copy.
11. **No axon model, no thresholds.** This stage produces the field only.

## Multi-lead FAST-style poster workflow

The macros now discover all Lead Designer leads in the active document, preserve
lead/contact identity, and support matched total-current scaling. The Potential
Visualizer includes editable waveform frequency (90 Hz default), phase width,
assumed gap, phase snapshots, and a self-contained animation/PNG export viewer.
See [poster workflow](../docs/poster_fast_workflow.md) for instructions and limits.

## SCS Job Manager

Open `SCS_JobManager.FCMacro` in FreeCAD to configure local mesh threads and
Elmer resources and monitor jobs across models. Updated modeling macros use
these shared execution profiles. The HPC profile (an SSH host plus sbatch
options) runs Generate Mesh and Solve / Plot Field on the cluster and fetches
the results back. See
[Job Manager workflow and limitations](../docs/job_manager.md).
