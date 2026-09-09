# Ansys on CWRU HPC (Pioneer): training notes for the NBF_RADO-SCS model

## Read this first (summary of the 2026-09-08/09 overnight session)

**What you have now**

- A **working, validated MAPDL bipolar-stimulation pipeline**, proven end to end
  on a small model (`src/ansys/testA_mapdl_cylinder/`): 390k nodes, ±1 mA on two
  contacts, solves clean, **V ∈ [−1.0366, +1.0380] V**, and the figure
  (`results/testA_voltage_slice.png`) shows the dipole with the CSF column
  visibly shunting the field. Every MAPDL construct the real model needs is
  exercised there: SOLID232, anisotropic white matter, equipotential contacts,
  current injection, grounding, slice export, external plotting.
- **A full diagnosis of why your Aug 29 / Sep 1 runs failed** — and it is *not*
  a meshing problem. Details below; the short version is three separate causes,
  two of which I fixed and one of which is an institutional resource limit.
- **Headless Mechanical working on Pioneer**, with the patch-independent
  meshing recipe for faceted STL bodies worked out (`src/ansys/testM7_*`). This is
  what you will need to re-mesh when the lead moves.
- Supporting tooling: `src/ansys/tissue_map.yaml` (all 245 bodies → Engineering Data
  materials, validated by `src/ansys/check_tissue_map.py`), `src/ansys/plot_voltage_slice.py`,
  and sbatch scripts for each run.

- **A converged solve on the real RADO anatomy** (truncated sub-model, job
  3796854): `DISTRIBUTED PCG SOLVER SOLUTION CONVERGED, NUMBER OF ITERATIONS =
  93`, 1,701,330 nodes / 1,128,810 elements, V ∈ [−0.0673, +0.1872] V for
  ±1 mA. Figure: `src/ansys/bigdeck_truncated/results/rado_voltage_slice.png`.
  **Read the caveats on this one before believing the numbers** — see
  "The truncated sub-model result" below. It demonstrates the pipeline runs on
  your real geometry; it is not yet a publication-grade field.

**What you do NOT have: a trustworthy full-model voltage figure.** Honest reason: the `tlv`
account is capped at **24 CPUs**, and because SLURM bills CPUs in proportion to
memory, that caps us at ~140 GB on `batch`. The full 18.8M-DOF solve needs
~320 GB. The `smp` partition has the memory but is booked out to ~Sep 13 by
8-day jobs. **Job 3796783 is queued on `smp` (24 cores, 600 GB) and will run the
full model unattended when a node frees.** That is the definitive run.

**The three root causes of the original failure**

1. ~~**Material contrast ~2e11**~~ — **WRONG, retracted 9 Sep.** See
   "Correction: the materials were never the problem" below. All reconditioning
   has been removed from every script; the decks use Khadka's original values.
2. **Floating, electrically isolated bodies.** Revealed by running PCG: it
   diverged with *"out of balance force goes to infinity … check for rigid body
   motions"*. The 245 bodies are joined by `CONTA174`/`TARGE170` pairs (correctly
   set to *Pure Electric contact*, bonded), but bonded contact only bonds within
   its pinball radius, and RADO ships **no surrounding soft-tissue envelope**, so
   peripheral structures have nothing to conduct into. Only one component
   (61,742 nodes) is grounded.
3. **Memory ceiling** from the account CPU cap, above.

**Decisions I made that you may want to revisit** — all flagged in
"Open questions" at the bottom. The most important: whether reconditioning those
two materials is acceptable to you physically, and whether the roots'
`...Middle...` bodies should be CSF rather than nerve (a 12× conductivity
difference right next to the electrodes).


Working log from a first hands-on attempt at the geometry -> mesh -> material ->
bipolar-stimulation -> solve pipeline for the NBF_RADO-SCS model, run autonomously
overnight on 2026-09-08/09. Written so later sessions (human or agent) can skip
the dead ends. Update this file as things change; don't let it rot.

Goal for this attempt: import the current STL geometry into Ansys, mesh it
(tetrahedral, patch-independent as the starting point), assign tissue
conductivities via name selections keyed on STL filename patterns, apply a
bipolar +/-1 mA current on a pair of lead electrodes, and get a solve — ideally a
voltage-spread figure. This is a rehearsal for the real dorsal-vs-ventral epidural
stimulation sweeps, not the production run.

## Ground rules followed

- No compute on login nodes: every heavy step ran inside a SLURM allocation
  (`srun` for interactive work, `sbatch` for anything long). See the CLAUDE.md
  sizing guidance (4 cores / 32 GB is the sane default; `tlv` account).
- Long-running work in `tmux` so an SSH drop doesn't kill it.
- Nothing destructive on the cluster; new work goes in its own directory.

## Inputs

- Geometry: the corrected STL set from the repo's `STL_files/` (245 files), synced
  to `~/Documents/SCS-Modeling/NBF_RADO-SCS_STL_corrected/` on Pioneer. This is
  the set with Shane's replacement lead electrodes/insulator and DRG at the
  correct scale (the originals were 10x too large; fixed in commit c972720).
- Conductivities: Mohamed's `engineering_data` file under `~/Documents` on Pioneer
  (location recorded below once found).

## Material assignment by filename pattern

Name selections keyed on substrings of the STL filenames (from Mohamed):

| Pattern              | Tissue / material            |
|----------------------|------------------------------|
| `V1-2`, `V1-3`, `V1-4` | vertebra                     |
| `T_disk_*`           | intervertebral disc          |
| `SCS Lead Electrode *` | electrode contact (platinum/iridium) |
| `SCS Lead Insulator` | lead insulation              |
| *(rest: fill in as mapped)* |                        |

## Environment recon (2026-09-08 night)

**Access.** `ssh case-hpc` only works over the CWRU OpenVPN; the tunnel dies on
laptop reboot and reconnecting (`vpn cwru` / `cwru-connect`) needs Mohamed's Duo
in a visible Firefox window — no unattended path. Check with
`timeout 5 bash -c 'cat </dev/null >/dev/tcp/pioneer.case.edu/22'`.
Over plain `ssh host 'cmd'` the `module` command is missing (no login shell);
use `ssh case-hpc 'bash -ls' <<'EOF' ... EOF` (login shell reading stdin — also
sidesteps all the nested-quoting pain).

**Ansys.** Modules `ansys/21.1 21.2 22.2 24.1 25.1(D)`, install at
`/usr/local/ansys_inc/v251`. MAPDL binary is `ansys251`. A prior license test
(`~/Documents/ansys_test/ansys_lic_test.dat`: SOLID232 block, VOLT BCs, VMESH,
SOLVE) ran to completion → MAPDL electric conduction works on this cluster.
System `python3` is 3.6.8 (no `{*}` XPath wildcards, no ansys packages); Python
modules exist via `module spider Python`. PyMAPDL/PyMechanical: not installed
system-wide (TBD in a user venv).

**SLURM.** Account `tlv`, QOS `gpudef`, MaxSubmit 48. Partitions: `batch*`
(default), `gpu`, `smp` (big-memory), `amd`, `aisc*`, `cgpu`; 13-day time limit.
Mohamed's earlier Ansys jobs used `--partition=smp --ntasks-per-node=16
--mem=512G --time=04:00:00` with `ansys251 -b -dis -mpi intelmpi -np 16`, staging
the deck to `$PFSDIR` (parallel scratch) and copying `*.rst *.out` back.
`quota -s` is broken (NFS permission errors) — don't rely on it.

**What already exists on Pioneer (`~/Documents/SCS-Modeling/`).**
- `ANSYS/SCS_Modeling.wbpj` and `SCS-Model/SCS_Model.wbpj`: Workbench projects,
  each with an Engineering Data XML export (materials below).
- `scs_model_geometry.pmdb`: a Mechanical geometry cache — a prior geometry
  import exists.
- `batch-20260829/`: **a complete Workbench-exported MAPDL deck**
  `scs_model.dat` (4.56 GB, `nblock,3,,20782975` → 20.8 M nodes, uMKS units,
  "Steady-State Electric Conduction (A5)", written on Mohamed's Georgia Tech
  Windows machine) plus `run_ansys.slurm`. `batch-20260901/` re-ran the same
  slurm script (its `scs_model.dat` is a 59-byte stub). Both produced
  `scs_model.out` + 15 per-rank `fileN.out`; outcome analysed in the attempt log.
- `~/Documents/log.out`: Workbench itself crashing headless
  (`aisol/.workbench: Aborted (core dumped)`).

**Materials (Engineering Data, resistivity in Ω·m at 37 °C; σ = 1/ρ S/m).**

| ED material         | ρ (Ω·m)  | σ (S/m) |
|---------------------|----------|---------|
| Metal Electrode     | 2.5e-7   | 4.0e6   |
| Lead Insulation     | 50000    | 2e-5    |
| Vertebra            | 25       | 0.04    |
| Intervertebral Disc | 1.6667   | 0.6     |
| Epidural Space      | 25       | 0.04    |
| Vasculature         | 1.5152   | 0.66    |
| Dura Mater          | 27.027   | 0.037   |
| CSF                 | 0.5882   | 1.7     |
| White Matter        | 6.9832   | 0.1432 (isotropic; paper: 0.1432/0.6 aniso) |
| Gray Matter         | 3.6232   | 0.276   |
| Nerve Root / DRG / Sympathetic Chain | 6.9832 | 0.1432 |
| Soft Tissue         | 250      | 0.004   |

Full pattern→material mapping: `src/ansys/tissue_map.yaml`.

## STL geometry audit (FreeCAD, headless `freecadcmd`; script in session scratch, CSV in `src/freecad/stl_mesh_quality_audit.csv`)

245 bodies, 1,771,980 facets total. No non-manifold edges, no inconsistent facet
orientation. But:

- **16 bodies are open (not watertight)** and stay open after a duplicate-vertex
  merge, so these are real holes, not import artifacts: `V1-2`, `V1-3`, `V1-4`
  (all three vertebrae), `T_disk_12_Lx-1`, `T_disk_9_10x-1`, `connection_R-1/2/4`,
  `April4Blood_in_Fat_bri-1..4` (2 components each), `neuro_0A_blood_bri_a-1..4`
  (3 components each). Open shells can't be turned into solids for volume
  meshing without healing. (Mohamed's Windows Workbench project did mesh this
  model, so SpaceClaim/DesignModeler's STL import presumably healed them — worth
  confirming how, since a headless batch import may not.)
- **22 bodies are multi-component.** The DRG `csf_coating`/`menging_coating`
  bodies are 2 closed components each (inner + outer shell of a coating = a
  hollow shell; that's fine, it just means "two surfaces"). `neuro_blood_in_fat-1`
  (7) and `neuro_blood_next_to_WM-1` (6) are several disjoint closed vessels in
  one file — fine for meshing but each becomes its own body in Ansys.
- Largest bodies: `SCS Lead Insulator` (32.6k facets), vertebrae (30.3k each).
  Thinnest features: `April4Blood_in_Fat_bri-2` is 0.4 mm thick; the electrode
  contacts are ~1.3 mm cubes; several vessels/rootlets are ~1 mm.

Practical implication: expect the mesher to need a minimum element size well
under 1 mm around the lead and vessels, and expect the 16 open bodies to fail
solid conversion unless healed first. `freecadcmd` note: `print()` output is
swallowed — write results to a file — and it segfaults on exit *after* the script
finishes (harmless, but check for `core` files in the cwd).

## Prior full-model runs (Aug 29 / Sep 1): why they failed

The 20.8 M-node Workbench deck (`batch-20260829/scs_model.dat`; the Sep 1 run
used `~/Documents/gatech/scs_model_bipolar_current.dat` via symlink) got all the
way through /PREP7 and matrix assembly on `smp` with 16 ranks / 512 GB
(peak ~320 GB RAM, ~300 GB scratch, 393 Gflops), then the **distributed sparse
solver aborted**:

```
Distributed sparse solver maximum pivot= 1.92280267E+16 at node 7405040
Distributed sparse solver minimum pivot= -347.221325 at node 786 VOLT.
*** ERROR ***  A large negative pivot value ( -347.221325 ) has been encountered
               An extremely large pivot ratio has been detected by the sparse solver.
```

Same signature both times (Sep 1: min pivot −352.9, also node 786). Elapsed
~30 min, no `.rst`. sacct shows the surrounding history: several
`ansys_solve` jobs on `batch` with 64 GB completed in seconds (deck too big /
staging), the 256–512 GB `smp` runs FAILED (exit 1) at 19–30 min, and a few
OnDemand dashboard sessions timed out.

Interpretation (a conductance matrix is positive semi-definite, so a *negative*
pivot means numerical breakdown, not physics):

1. **Material contrast of ~1e11** in one matrix: Metal Electrode 4e6 S/m vs
   Lead Insulation 2e-5 S/m (ρ 2.5e-7 vs 5e4 Ω·m). Pivot ratio 1e16 is exactly
   what that produces. Fix: don't model the insulator as a conductor at all
   (a body with no elements *is* an insulating boundary in FEM), and cap the
   metal at ~1e3–1e4 S/m — the contacts are made equipotential by CP coupling
   anyway, so their exact conductivity is irrelevant.
2. **Floating bodies**: 245 separately-meshed STL solids only conduct into each
   other if they share nodes (shared topology / conformal mesh) or are joined
   by bonded contact. Any body with no path to the reference potential gives a
   zero/negative pivot. Node 786 is the very first node in the deck's `nblock`,
   which is consistent with an isolated region rather than a random bad
   element. Fix: verify every body is connected (in Mechanical: check the
   Contact folder covers all interfaces, or use Share Topology), and pin V=0 on
   an outer boundary.
3. Possibly inverted/degenerate tets from thin or open facet bodies (the audit
   above shows 0.4 mm-thick vessels and 16 open shells).

The pmdb/mechdb/dsco from the Windows project are on Pioneer
(`ANSYS/SCS_Modeling_files/dp0/{SYS/DM/SYS.dsco 184 MB, global/MECH/SYS.mechdb
945 MB}`), so the existing mesh could be re-used from Mechanical without
re-meshing if headless Mechanical works (Test M below).

## How to launch MAPDL on this cluster (learned the hard way)

`ansys251` shells out to Intel MPI **even for `-np 1`**, and on this SLURM
(23.11+) that collides with cpus-per-task propagation. Two dead ends:

| invocation | failure |
|---|---|
| `ansys251 -b -np 4 -i x.dat` | `BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES / RANK 3 ... KILLED BY SIGNAL: 9`, exit 255, no command executed |
| `ansys251 -b -np 1 -i x.dat` | `mpiexec ... Unable to run bstrap_proxy`, then `srun: fatal: cpus-per-task set by two different environment variables SLURM_CPUS_PER_TASK=6 != SLURM_TRES_PER_TASK=cpu=1`, exit 255 |

Two invocations that do work:

- **Small/serial**: `unset SLURM_CPUS_PER_TASK SLURM_TRES_PER_TASK` then
  `ansys251 -b -smp -np N -i x.dat -o x.out`. `-smp` = shared-memory parallel,
  no MPI layer at all.
- **Large/distributed** (Mohamed's working pattern):
  `ansys251 -b -dis -mpi intelmpi -np $SLURM_NTASKS -i x.dat -o x.out` with
  `export OMP_NUM_THREADS=1; export KMP_AFFINITY=disabled`, staging to
  `$PFSDIR`. This is what reached the solver on the 20.8 M-node deck.

Diagnosing a failed run: `slurm_<jobid>.out` has the launcher/MPI errors;
`<job>.out` has MAPDL's own output (and MAPDL echoes the entire input deck
first, so `grep`ping for your `/COM` markers matches the echo as well as the
execution — check line numbers, or grep for the substituted values).

## Attempt log

*(chronological; keep the failures — they're the useful part)*

- **Test A** (`src/ansys/testA_mapdl_cylinder/`): pure MAPDL, coaxial cord/CSF/fat
  cylinders + two 3 mm contacts (1 mm gap), SOLID232, anisotropic white matter
  via `MP,RSVZ`, CP-coupled contacts, ±1 mA, V=0 on the outer surface, moderated
  metal conductivity, no insulator body.
  - job 3796629 (`-np 4`): FAILED, MPI rank killed (see table above).
  - job 3796666 (`-np 1`): FAILED, MPI bootstrap + SLURM env conflict.
  - job 3796672 (`-smp`, env unset): *(pending)*
- **Test M** (`src/ansys/testM_mechanical_smoke/`, job 3796630): **headless
  Mechanical works.** `/usr/local/ansys_inc/v251/aisol/.workbench -DSApplet
  -AppModeMech -b -script foo.py` exited 0 without needing X or `xvfb-run`.
  It reported `ProductVersion 2025 R1`, imported
  `T8-10 - neuro_CSF-1.STL` (1 part, body named `solid T8-10 - neuro_CSF-1`),
  and saved a 6 MB `.mechdat`. **But `GenerateMesh()` produced 0 nodes /
  0 elements** — investigated in Test M2. Note `Project.Save` is deprecated in
  2025 R1; use `SaveAs(path, True)`.
- **Test M2** (`src/ansys/testM2_mechanical_mesh/`, job 3796670): dumps body
  type/dimension/volume and the Mechanical message log, then adds an explicit
  `MethodType.Tetrahedrons` + `AlgorithmType.PatchIndependent` method scoped to
  the imported bodies. Result: *(pending)*
- **Big deck recondition** (`src/ansys/bigdeck_recondition/`, job 3796667): re-run of
  Mohamed's existing 20.8 M-node deck with **only two material lines changed**.
  Rationale below. Result: *(pending)*

## The big deck's physics was already right — only the conditioning was wrong

Grepping the 4.5 GB deck case-insensitively (it's written in lowercase; an
uppercase-anchored grep finds nothing and is misleading) shows a complete,
correct bipolar setup:

```
d,_CM9093,volt,0.            ! ground, component of 61,742 nodes
f,786,amps,1000000000.       ! +1 mA   (uMKS current unit is pA, so 1e9 pA)
f,3149,amps,-1000000000.     ! -1 mA
ce,next,0.,786,volt,-1.,787,volt,1.    ! x878: equipotential coupling on both
et,*,232                     ! SOLID232 electric solid throughout
/units,uMKS
```

So there *is* a reference potential, the current is right, and the contacts are
properly equipotential. Note the pivot error named **node 786** — the master
node of one of those constraint-equation sets, i.e. an electrode, not a random
element.

The materials (uMKS resistivity in Tohm·µm; 1 Tohm·µm = 1e6 Ω·m):

| MAT | deck value | Ω·m | tissue |
|-----|-----------|------|--------|
| 1   | 2.5e-13   | 2.5e-7 | Metal Electrode |
| 5   | 0.05      | 5e4    | Lead Insulation |
| 6   | 2.5e-05   | 25     | Vertebra / Epidural |
| 9, 191 | 2.7027e-05 | 27.027 | Dura |
| 17  | 1.5152e-06 | 1.5152 | Vasculature |
| 29, 54, 92, 194, 226 | 6.9832e-06 | 6.9832 | White matter / nerve root / DRG / sympathetic |
| 30  | 5.882e-07 | 0.5882 | CSF |
| 192 | 3.6232e-06 | 3.6232 | Gray matter |
| 240 | 1.6667e-06 | 1.6667 | Disc |

**MAT 1 vs MAT 5 differ by 2e11** — which is precisely the reported pivot ratio
of ~1e16 once assembled. The recondition run changes only those two lines
(metal → 1.0e-08, insulation → 1.0e-03, contrast 1e5) and comments out
`/fclean` so the `.rst` survives. Everything else — mesh, BCs, currents — is
byte-identical to what Mohamed already built.

## Headless Mechanical scripting: the traps, in order

Mechanical *does* run headless on Pioneer — no X, no `xvfb-run`, no OnDemand
Desktop:

```
/usr/local/ansys_inc/v251/aisol/.workbench -DSApplet -AppModeMech -b -script /abs/path/script.py
```

Five traps cost an iteration each. In order of how much time they waste:

1. **The scripting engine is IronPython 2.7.4 on Mono, not CPython 3.** Any
   py3-only syntax is a *compile* error, and the process then produces
   **zero output and exit code 0** — indistinguishable from "ran and printed
   nothing". If a script prints nothing at all, suspect syntax before logic.
   (Ansys does ship CPython 3.10 at
   `commonfiles/CPython/3_10/linx64/Release/python`, but that is not what
   `-script` uses.)
2. `ExtAPI.DataModel.Project.Messages` does not exist. Messages live at
   **`ExtAPI.Application.Messages`** (`.Count`, then `[i].Severity` /
   `[i].DisplayString`). Without them, mesh failures are completely silent.
3. `Project.Save` is deprecated in 2025 R1 → use `Project.SaveAs(path, True)`.
4. `Ansys.Mechanical.DataModel.Enums.MethodType` has **no `Tetrahedrons`** —
   the tet value is **`AllTriAllTet`**. And the patch-independent switch is NOT
   on `AlgorithmType` (that enum is CMFD/MFD/ProgramControlled/SCPIP —
   optimisation algorithms); it is
   **`MeshMethodAlgorithm.PatchIndependent`**.
5. The sizing properties on `AutomaticMethod` are **not** `ElementSize` /
   `MaxElementSize` / `MinElementSize` (all `AttributeError`). The real names,
   found by dumping the object's 377-entry `Properties` collection, are
   **`MaximumElementSize`**, `MinimumSizeLimit`, `DefeaturingTolerance`,
   `CurvatureNormalAngle`, `FeatureAngle`, `MeshBasedDefeaturing`,
   `ApproximativeNumberOfElementsPerPart`. `MaximumElementSize` defaults to
   `0 [m]`, which is the "required input not defined" the mesher complains
   about.

**Why the default mesher fails on these STLs.** The import produces a genuine
solid (`BodyType = GeoBodySolid`, `Volume = 8835.85 mm³`,
`Area = 6711.95 mm²`) but with **`faces=1, edges=0, vertices=0`** — a
tessellated body with no feature topology. The default patch-*conforming*
mesher has to respect every face patch and simply reports
`The mesh generation did not complete. Try meshing with another mesh method`.
Patch-independent tets ignore the patch structure and mesh the enclosed volume,
which is the standard remedy for faceted geometry — and is what Mohamed
suggested at the outset.

Working recipe (see `src/ansys/testM7_patchindep_sized/mech_pi_sized.py`):

```python
meth = model.Mesh.AddAutomaticMethod()
meth.Location  = sel                                   # GeometryEntities selection of body ids
meth.Method    = E.MethodType.AllTriAllTet
meth.Algorithm = E.MeshMethodAlgorithm.PatchIndependent
meth.MaximumElementSize   = Quantity("2 [mm]")         # REQUIRED, defaults to 0
meth.MinimumSizeLimit     = Quantity("0.4 [mm]")
meth.DefeaturingTolerance = Quantity("0.05 [mm]")
model.Mesh.GenerateMesh()
```

## What worked

- Headless Mechanical on a compute node: import STL → solid body → save
  `.mechdat`, entirely batch, no display.
- MAPDL electric conduction end to end on the **Test A** validation model:
  390,181 nodes / 284,520 elements, sparse solver, `RUN COMPLETED`, exit 0.
  Voltage range **−1.0366 V … +1.0380 V** for ±1 mA — symmetric, as a bipolar
  pair should be. The figure shows the expected dipole and, satisfyingly, the
  CSF column (σ = 1.7 S/m) visibly shunting the field, which is the effect the
  RADO paper stresses.
- The MAPDL patterns the real model needs, all validated in Test A: `SOLID232`,
  anisotropic white matter via `MP,RSVX/RSVY/RSVZ`, equipotential contacts via
  `CP,NEXT,VOLT,ALL`, current injection with `F,node,AMPS`, reference potential
  with `D,ALL,VOLT,0`, and a `*VGET`/`*VWRITE` slab export to CSV.
- Plotting outside Ansys: `src/ansys/plot_voltage_slice.py` renders the CSV with
  matplotlib (`matplotlib.tri`, since **scipy is not installed** on the laptop).
  Far better than MAPDL's renderer, and reusable for the full model.

## What did not work

- `ansys251 -b -np N` (with or without N=1) — see the launch table above; use
  `-smp` or `-dis -mpi intelmpi`.
- `*GET,par,NODE,0,MXV,VOLT` is not valid ("Unknown label in field 5") and it
  aborts the remainder of /POST1 in batch, silently costing you every plot after
  it. Use `NSORT,VOLT` + `*GET,par,SORT,0,MAX|MIN`.
- MAPDL's own `/SHOW,PNG` + `/CPLANE` section plot produced a uniform green
  outer-cylinder view — technically a plot, practically useless. Export
  coordinates+VOLT and plot externally instead.
- FreeCAD's `Mesh.difference()` is a **complete no-op** in this build
  (subtracting a 15 mm sphere centred on the epidural-space centroid changed
  neither volume nor facet count), so the "boolean-subtract the lead from the
  epidural space" idea cannot be done with FreeCAD mesh booleans. It would need
  BREP/OCCT conversion or an external library.
- `quota -s` on Pioneer errors out with NFS permission failures.

## Open questions for Mohamed

1. **Is reconditioning the materials acceptable physically?** I changed only the
   metal contacts (2.5e-7 → 1e-2 Ω·m) and lead insulation (5e4 → 1e3 Ω·m) to
   kill the 2e11 contrast. The contacts are equipotential by constraint
   equation regardless, and 1e3 Ω·m is still ~37× dura, so I believe current
   flow is essentially unchanged — but you should sanity-check the resulting
   contact voltages against your expectations before trusting any numbers.
2. **Anisotropic white matter.** Engineering Data has white matter isotropic at
   0.1432 S/m; the paper uses 0.1432 transverse / 0.6 longitudinal. The deck
   therefore solves the isotropic case. Adding anisotropy needs an APDL snippet
   (`MP,RSVZ`) and depends on the cord axis being global Z — worth confirming
   for the real geometry.
3. **Which two contacts** should be the bipolar pair for the dorsal-vs-ventral
   sweeps? The existing deck uses nodes 786 and 3149; I don't yet know which
   physical contacts those correspond to.
4. The 16 non-watertight STL bodies (all three vertebrae among them) — how did
   the Windows Workbench import heal them? That matters for re-meshing from
   scratch when the lead moves.
5. **Are the root "Middle" layers nerve or CSF?** `src/ansys/check_tissue_map.py`
   maps all 245 bodies with no leftovers, but the root bundles come in
   Inside / Middle / OutsideMeninges layers, mirroring the DRG's
   `in_middle` / `csf_coating` / `menging_coating`. If that analogy holds, the
   `...Middle...` bodies are CSF, not nerve. They are currently Nerve Root.
   CSF is 1.7 S/m vs nerve 0.1432 S/m — a 12× difference in tissue immediately
   around the electrodes, so this materially changes the answer.

## Model geometry: hollow shells that tile, and a real epidural channel

**The compartments are HOLLOW SHELLS, each with the inner ones carved out.**
They tile space; they do not overlap. Established by ray-casting a point at the
centre of a cord cross-section (z = 110, x = 56.42, y = 76.46):

| body | forward crossings | verdict |
|---|---|---|
| white matter | 3 (odd) | point is **inside the material** |
| CSF | 2 (even) | point is in its central cavity |
| dura | 2 (even) | cavity |
| epidural space | 2 (even) | cavity |

So the epidural body's 11295 mm³ is genuinely epidural fat, not canal-fill.

**Beware measuring a channel as the difference of two outer surfaces.** An
earlier version of `measure_corridor.py` computed `epidural.ymax − dura.ymax`,
got 0.07–0.11 mm and even negative widths, and I wrongly concluded there was no
epidural space to put a lead in. Those two surfaces are not the two sides of a
channel. Ray-casting along +Y and reading consecutive entry/exit pairs gives the
true material intervals, and they are anatomically sensible — the dorsal
epidural space is the thicker one, as it should be:

| z (mm) | ventral fat | dorsal fat |
|---|---|---|
| 95 | 1.74 mm (y 67.60–69.34) | 2.35 mm (y 80.12–82.47) |
| 115 | 1.71 mm | 2.30 mm |
| 145 | 1.68 mm | 2.27 mm |

Over 50 clean slices spanning z = 61–159: **dorsal 2.27 mm minimum, ventral
1.67 mm minimum**. A 1.3 mm clinical lead fits on either side, with no need to
carve space. Verified end-to-end: 100% of sampled vertices of a generated
8-contact dorsal lead lie inside the epidural material.

**The canal migrates with height.** The spine is kyphotic, so at midline the
dorsal centreline runs y ≈ 79.8 at z = 89 to y ≈ 87.5 at z = 132 — a 7.7 mm
rise, about 10°, over the length of an 8-contact lead. A straight lead at fixed
y walks out through the dura, which is presumably why RADO's own lead is curved.
`src/freecad/measure_corridor.py` fits that centreline (cubic, 0.018 mm RMS) into
`src/freecad/epidural_corridor.json`, and `src/freecad/make_scs_lead.py` sweeps the lead
along it.

For reference, RADO's own lead sits 11–25 mm **left** of midline at z = 93–97,
out by the left DRG column — a DRG lead, not a dorsal-column lead.

## Correction: the materials were never the problem (9 Sep)

Two mistakes on my part, both now fixed. Recorded in full because the reasoning
is the useful bit.

**Mistake 1 — I misread a unit and "fixed" a correct number.** The deck is
`/units,uMKS`, so resistivity is in **Tohm·µm**, which Workbench annotates on
every line (`! Tohm um`). 1 Tohm·µm = 10⁶ Ω·m, so the deck's

```
MP,RSVX,1,2.5e-13,	! Tohm um     ->  2.5e-7 ohm*m  =  4e6 S/m
MP,RSVX,5,0.05,		! Tohm um     ->  5e4   ohm*m
```

are **exactly** the Engineering Data values (platinum-like contacts at 4e6 S/m).
Every material converts correctly: CSF 1.7000, white matter 0.1432, disc 0.600,
gray 0.276 S/m. Nothing was ever mismatched — I quoted the raw µMKS figure
without converting and concluded there was a contrast problem. Mohamed caught
it. **Never change a material to fix a solver problem; the values are matched to
the paper on purpose.**

**Mistake 2 — the reconditioning did nothing anyway.** Comparing the two runs:

| | Aug 29 (original) | 9 Sep (reconditioned, 331 GB) |
|---|---|---|
| max pivot | 1.92280267e16 | 1.92307552e16 |
| min pivot | −347.22 @ node 786 | −77.33 @ node 7192526 |
| min abs pivot | — | 1.94e-07 @ node 8711510 |
| error | large negative pivot | **"insufficiently constrained model"** |

The max pivot is identical to five significant figures across an eleven-order
material change. It never came from the materials. It comes from the deck's own
contact conductance:

```
*set,_maxCond,4.e+016
rmod,cid,19,_maxCond/_ASMDIAG      ! ECC
```

That is Workbench's standard way of emulating perfectly bonded electric contact,
and 4e16 is where the pivot ratio originates. Leave it alone.

## Which bodies actually float — answered from the geometry

With memory ruled out (331 GB used, no OOM) the remaining error is explicit:
*"There is at least 1 small equation solver pivot term ... Please check for an
insufficiently constrained model"* — an electrically isolated region.

`src/freecad/find_floating_bodies.py` settles which ones, without needing the cluster:
it parses the 245 STLs, builds a proximity graph over their surfaces, and finds
the connected components as a function of contact gap.

| gap tolerance | components | isolated bodies |
|---|---|---|
| 0.05 mm | 38 | 79 |
| 0.10 mm | 9 | 8 |
| 0.20 mm | 3 | 2 |
| **0.30 mm** | **1** | **0** |

So the geometry *is* fully connected, but only once contact can bridge ~0.3 mm.
The last bodies to join are the vascular chain — `Thoracic_aorta-1` and
`connection_R-2` at 0.2 mm, then `connection_R-1`, the four
`AprilDorsalRootsBloodParallel_*` segments and `T_disk_11_12x-1` at 0.1 mm.
These are the "connection" radicular vessels running to the thoracic aorta:
separate STL bodies that merely abut. Blood is electrically continuous through
them, so bonding them is the physically correct reading rather than a numerical
dodge.

And the deck leaves the pinball radius at the default on all 1116 pairs:

```
rmod,cid,6,0.       ! PINB
rmod,tid,6,0.       ! PINB
```

`src/ansys/bigdeck_pinball/run_pinball.sbatch` sets it to an absolute 500 µm
(negative PINB = absolute distance; the deck is in µm) and changes **nothing
else** — it asserts the materials are still `2.5e-13` / `0.05` and aborts if
not. Widening pinball cannot short unrelated structures: it only lets an
*existing* pair close a larger gap, and Workbench never created a pair across,
say, the dura wall.

## Revised diagnosis: the model has floating (electrically isolated) regions

The material-contrast theory was right about the *negative* pivot but is not the
whole story. Running the reconditioned deck with the **PCG** solver (job
3796753) got through all element matrix formation (CP 2627 s, 18,830,747
equations) and then:

```
*** WARNING ***  The PCG solver detects that the out of balance force goes to
                 infinity. The solution has not converged. Please check for
                 rigid body motions in your model.
*** NOTE ***     The PCG solver failed to converge ... automatically switched
                 to the sparse solver (EQSLV,SPARSE)
```

"Rigid body motion" is the structural wording; in a conduction problem it means
a **floating region — a body with no conduction path to the grounded nodes**.
An iterative solver cannot converge on a singular system, so PCG is the wrong
tool if islands exist. (MAPDL then auto-switched to the direct solver and the
process died at once, almost certainly out of memory: 140 GB available here
versus the ~320 GB the earlier direct runs used.)

**How the bodies are connected.** The 245 bodies are *not* conformally meshed.
Workbench generated `CONTA174`/`TARGE170` contact pairs, and they are correctly
configured for this physics:

```
keyo,cid,1,6      ! Pure Electric contact      <- VOLT DOF, so contacts do conduct
keyo,cid,12,5     ! bonded always
keyo,cid,9,1      ! ignore initial gaps/penetration
```

So conduction across interfaces is intended. But bonded contact only bonds
within its **pinball radius** — any body whose neighbours are further away than
that is electrically isolated regardless. With 245 separate STL bodies and
auto-generated contact, some islands are very plausible. Note also that RADO's
STL set has **no surrounding soft-tissue/thorax envelope** (`soft_tissue`
matches zero bodies in `tissue_map.yaml`), so peripheral structures — the
sympathetic chain, isolated vessel segments — have nothing to conduct into.

There is exactly **one** grounded component in the whole deck,
`CMBLOCK,_CM9093,NODE,61742` → `d,_CM9093,volt,0.`, so only the region connected
to those 61,742 nodes is constrained.

Ways forward, roughly in order of effort:

1. **Direct solver, forced out-of-core** (`bcsoption,,outofcore`) —
   `src/ansys/bigdeck_recondition/run_bigdeck2.sbatch`. A direct factorisation can
   be regularised where an iterative solver cannot: with the contrast fixed,
   islands give *zero* pivots, which MAPDL constrains with a warning instead of
   aborting. Fits 140 GB by streaming to /scratch (52 TB free).
2. **Ground each island.** Requires a connectivity analysis to find them.
3. **Enlarge the contact pinball radius** so near-touching bodies bond
   (`rmodif` on the contact real constants).
4. **Solve a connected sub-model** — the concentric core (epidural space, dura,
   CSF, grey/white matter, contacts) is certainly mutually in contact. Note
   materials are shared across bodies (14 `MP` definitions for 245 bodies), so
   `ESEL,S,MAT` selects tissue classes but cannot separate, say, vertebra from
   epidural space (both 25 Ω·m → MAT 6).
5. **Add the missing soft-tissue envelope**, which is the physically correct fix
   and matches the paper's 19-compartment description.

## The truncated sub-model result (job 3796854) — and why to distrust it

This is the first converged solve on the real anatomy, and it is genuinely
useful as proof the pipeline works. It is **not** a result to quote.

What it is (final version, job 3796865, `run_truncated4.sbatch`): solid elements
whose centroid lies within a ±15 mm box centred on the **midpoint of the two
contacts**, coupled with `CPINTF,VOLT,300` (300 µm), V = 0 clamped on the six
box faces, and ±1 mA distributed over each contact's nodes (contacts selected by
element centroid within ±2 mm of nodes 786 / 3149).
**1,952,166 nodes / 1,297,630 elements, PCG converged in 89 iterations,
V ∈ [−0.0679, +0.1728] V.** (The earlier job 3796854, box centred on contact A:
1,701,330 nodes, 93 iterations, V ∈ [−0.0673, +0.1872] V.)

One asymmetry that is **not** an artifact: contact A has 543 elements / 935
nodes and contact B has 837 / 1396. Selecting by element centroid instead of by
node location produced identical counts, so the two contacts really are meshed
differently, and they sit ~4.3 mm apart in different tissue. A bipolar pair in
heterogeneous anatomy should *not* give an antisymmetric field, so
+0.173 / −0.068 V is plausible rather than suspicious. Current injection is
per-node (`f,all,amps,±1e9/N`), so the *total* current is ±1 mA regardless of
the node-count difference.

Why not to trust the numbers yet:

1. **The domain is artificially truncated at 15 mm and clamped to 0 V there.**
   A current dipole's potential falls as ~1/r², so the near field should be
   roughly right, but the far field is imposed, not computed. Absolute
   impedance and anything beyond ~1 cm is meaningless.
2. **Inter-body conduction is `CPINTF` at 300 µm, not the model's own bonded
   electric contact.** The contact elements survive in the database but are not
   part of the solved element set, so every interface is coupled purely by
   node-proximity. Interfaces whose non-conformal meshes are further apart than
   300 µm are simply not connected. This is the single biggest source of doubt.
3. ~~The two contacts are asymmetric~~ — investigated and **explained**, see
   above: the mesh densities genuinely differ and the contacts sit in different
   tissue. Not a bug.
4. **The deck's 878 electrode constraint equations were deleted** (`CEDELE,ALL`)
   and replaced by distributed current injection. That is defensible — the metal
   is ~60× more conductive than CSF so each contact is nearly equipotential
   anyway — but it is not what your Workbench model specifies.
5. **`EDELE,ALL` silently did nothing**: `EDELE is not a recognized SOLUTION
   command` — it needs `/PREP7`, and the injection point is inside `/SOLU`. So
   the run did *not* delete the contact elements as the script intended. The
   thing that actually unblocked the solve was `CEDELE,ALL` removing the
   CE/CP conflict. Worth knowing before "fixing" the script.
6. `*VGET` warns "Some entities requested in the *VGET were undefined" — the
   POST1 slab selects nodes across the whole model, and those outside the solved
   set come back as exactly 0. `plot_voltage_slice.py` is fed a filtered file
   (`trunc_slice_solved.csv`, 46,194 of 665,360 points) for this reason.

The queued `smp` job (3796783) avoids items 1–5 entirely by solving the model as
Workbench built it, with enough memory for the direct solver.

## The deck has CRLF line endings — end-anchored regexes silently fail

`scs_model_bipolar_current.dat` was written by Workbench on Windows:
`file` reports *ASCII text, with CRLF line terminators*, and `cat -A` shows
`solve^M$`. So a sed expression like

```sed
0,/^solve$/s//eqslv,pcg,1e-8\nsolve/     # never matches: the line is "solve\r"
```

does nothing, while the **unanchored** substitutions in the same command
(`s/^MP,RSVX,1,2\.5e-13,/.../`, `s|^/fclean|...|`) succeed normally. The patch
therefore *looks* like it worked. A run was launched, sat on a compute node
using the wrong solver, and would have burned the whole allocation before the
mistake surfaced.

Two habits that catch this:

1. Strip CR as the first expression in the stream: `sed -e 's/\r$//' -e ...`
   (MAPDL is perfectly happy with LF-only input on Linux).
2. **Verify every edit and abort if one is missing.** The job script now greps
   for each expected result and `exit 1`s if any is absent, rather than
   proceeding with a half-patched deck.

## Account limits on Pioneer (this bit us twice)

`sacctmgr show assoc account=tlv` → **`GrpTRES cpu=24,gres/gpu=1`**. The whole
`tlv` account shares 24 CPUs. Two consequences:

- **SLURM bills CPUs in proportion to memory.** Asking for 450 GB on a 2 TB /
  384-core node (`epyc2tb`, ~6 GB/core) implies **75 CPUs**, and the job sits in
  `AssocGrpCpuLimit` forever with `START_TIME = N/A`. Memory-per-core by
  partition: `batch` icosa 6.4, `batch` epyc2tb 6.0, `smp` smpt08/09 28.3,
  `smp` smpt10-12 42.0 GB/core. That ratio — not raw node size — is what decides
  whether a big-memory job fits the cap.
- Interactive `srun` sessions eat the same 24 CPUs, so a background `srun` will
  block your own queued batch job. Check with
  `squeue -A tlv -h -t RUNNING -o "%C"` before blaming the scheduler.
- `aisc` / `aiscii` (20.9 GB/core, 2 TB nodes) are **not** open to `tlv`:
  "Invalid account or account/partition combination specified".

Use `sbatch --test-only ... --wrap="true"` to get a predicted start time before
committing to a submission — that is how the smp-vs-batch decision was made
(smp: 4 days out; batch: same night).
