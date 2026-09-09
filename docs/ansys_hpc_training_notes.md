# Ansys on CWRU HPC (Pioneer): training notes for the NBF_RADO-SCS model

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

Full pattern→material mapping: `ansys/tissue_map.yaml`.

## STL geometry audit (FreeCAD, headless `freecadcmd`; script in session scratch, CSV in `ansys/stl_mesh_quality_audit.csv`)

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

- **Test A** (`ansys/testA_mapdl_cylinder/`): pure MAPDL, coaxial cord/CSF/fat
  cylinders + two 3 mm contacts (1 mm gap), SOLID232, anisotropic white matter
  via `MP,RSVZ`, CP-coupled contacts, ±1 mA, V=0 on the outer surface, moderated
  metal conductivity, no insulator body.
  - job 3796629 (`-np 4`): FAILED, MPI rank killed (see table above).
  - job 3796666 (`-np 1`): FAILED, MPI bootstrap + SLURM env conflict.
  - job 3796672 (`-smp`, env unset): *(pending)*
- **Test M** (`ansys/testM_mechanical_smoke/`, job 3796630): **headless
  Mechanical works.** `/usr/local/ansys_inc/v251/aisol/.workbench -DSApplet
  -AppModeMech -b -script foo.py` exited 0 without needing X or `xvfb-run`.
  It reported `ProductVersion 2025 R1`, imported
  `T8-10 - neuro_CSF-1.STL` (1 part, body named `solid T8-10 - neuro_CSF-1`),
  and saved a 6 MB `.mechdat`. **But `GenerateMesh()` produced 0 nodes /
  0 elements** — investigated in Test M2. Note `Project.Save` is deprecated in
  2025 R1; use `SaveAs(path, True)`.
- **Test M2** (`ansys/testM2_mechanical_mesh/`, job 3796670): dumps body
  type/dimension/volume and the Mechanical message log, then adds an explicit
  `MethodType.Tetrahedrons` + `AlgorithmType.PatchIndependent` method scoped to
  the imported bodies. Result: *(pending)*
- **Big deck recondition** (`ansys/bigdeck_recondition/`, job 3796667): re-run of
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

## What worked

## What did not work

## Open questions for Mohamed
