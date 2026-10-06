# Agent coordination — `cad-dev`

Working file for the CAD rebuild of the RADO-SCS geometry. **Mohamed relays
between agents; neither agent sees the other's session.** Every task packet
below is therefore written to be executed cold, with no prior context.

The previous version of this file was a 774-line chronological log of the
September poster push. It is archived at `AGENTS_log_archive_2026-09.md`
(gitignored) and nothing in it is needed to work from this file.

---

## Roles

| | owns | does not touch |
|---|---|---|
| **Claude** (coordinator) | task definition, acceptance criteria, review of returned work, `fem/**`, `src/neuron/**` | CAD bodies while a task is out with Codex |
| **Codex / astra** | the CAD rebuild on `cad-dev`: `src/freecad/geometries/**` (new), geometry pilots | `fem/**`, `src/neuron/**`, `src/freecad/lead_designer.FCMacro`, `src/freecad/build_lead_config.py` |
| **Mohamed** | relay, priorities, anatomical sign-off | — |

## Hard constraints — apply to every task

1. **Work on branch `cad-dev`.** It exists. Do not commit to `main`.
2. **Do not modify the existing FEM or NEURON pipelines.** `fem/` and
   `src/neuron/` are a working, validated route; the CAD rebuild must not
   perturb them. New CAD work goes in a new directory.
3. **Never save a `.FCStd` the user has open**, and never run FreeCAD document
   edits under `freecadcmd` (it destroys per-face colours). Mohamed usually has
   `NBF_RADO-SCS_poster_dorsal.FCStd` open.
4. **Do not delete or overwrite anything under `STL_files/`.** It is the only
   copy of the source geometry in this repo.
5. Report numbers you measured, not numbers you expected. If a task can't be
   completed as specified, say which part and why rather than substituting a
   different deliverable.

## Context Codex will not otherwise have

- Repo root `/home/mohamed/Projects/SCS-Modeling`. FreeCAD 26.3.dev.
- `STL_files/` holds **245 STLs** named `T8-10 - <structure>-<n>.STL`. They are
  the RADO-SCS 3.0 release geometry (Khadka et al. 2020, J. Neural Eng. 17:026033).
- `src/ansys/tissue_map.yaml` maps those files into **11 tissue classes**; the
  mapping is read by `fem/scripts/config.py` as `TISSUE_BODIES`.
- **The central known defect of this geometry:** the STL compartments are *not*
  tessellated compatibly at shared interfaces — dura and epidural share 4
  vertices out of ~1780. They are individually watertight and manifold;
  compatibility is the problem. No conformal volume mesh can be built by
  stitching them, which is why the FEM route sidesteps them entirely with a box
  mesh plus point-in-shell tet classification. **This is the problem the CAD
  rebuild exists to solve.** Do not spend time rediscovering it.
- The model spans three vertebral bodies, labelled T10/T11/T12 in this repo
  (the RADO release calls the same anatomy T9–T11; the filenames say T8-10 —
  a known three-way naming inconsistency, not yet reconciled).

## Acceptance criterion for any reverse-engineered body

Derived from the pipeline, not chosen arbitrarily. The solved poster meshes used
`H_MIN = 0.2 mm`, so geometric error must sit well inside one element:

- **RMS deviation from the source mesh ≤ 0.05 mm**
- **Max (Hausdorff) deviation ≤ 0.20 mm**
- Body is watertight, manifold, and reports a positive volume within **1%** of
  the source mesh's enclosed volume.

A body that cannot meet this is a finding, not a failure — report it with the
numbers and move on to the next.

---

## Task board

| # | task | owner | state |
|---|---|---|---|
| 1 | STL triage inventory | Codex | **OPEN** |
| 2 | Reverse-engineering pilot, 3 bodies | Codex | blocked on #1 |
| 3 | Rostrocaudal tiling feasibility | Claude | not started |

---

### TASK 1 — STL triage inventory  *(Codex, open)*

**Why.** Before committing to rebuilding 245 bodies, we need to know what kind
of shapes they are. A revolve-like vertebral body and a freeform rootlet need
completely different treatment, and some may not be worth rebuilding at all.

**Do.** Produce `src/freecad/geometries/stl_inventory.csv` plus a short
`src/freecad/geometries/INVENTORY.md`, covering every file in `STL_files/`:

| column | meaning |
|---|---|
| `filename` | as on disk |
| `tissue_class` | from `src/ansys/tissue_map.yaml` |
| `n_facets`, `n_vertices` | raw counts |
| `watertight`, `manifold` | booleans |
| `volume_mm3`, `bbox_x/y/z_mm` | enclosed volume and extents |
| `genus` | 0 = topological ball, >0 = has handles |
| `planar_fraction` | fraction of facet area lying on detectable planes |
| `cylindrical_fraction` | fraction fitting a cylinder/revolve within 0.05 mm |
| `shape_hint` | one of `prismatic`, `revolve`, `swept`, `freeform` |

Then in `INVENTORY.md`, 1–2 pages: how many bodies fall in each `shape_hint`,
which tissue classes are dominated by which hint, the five largest and five
smallest bodies by facet count, and **your recommendation for which 3 bodies
Task 2 should pilot** — ideally one clearly prismatic/revolve (a vertebra), one
smooth closed surface (the cord or dura), and one hard case (a rootlet or
vessel).

**Constraints.** Read-only with respect to `STL_files/`. Pure Python is fine
(`numpy-stl`, `trimesh`) — this task does not need FreeCAD. Do not attempt any
reconstruction yet.

**Done when** both files exist on `cad-dev`, the CSV has 245 rows, and the
recommendation names three specific files with a reason each.

---

### TASK 2 — Reverse-engineering pilot, 3 bodies  *(Codex, blocked on #1)*

Scoped deliberately to three bodies. Do not start on the other 242 until the
pilot's numbers are reviewed.

**Do.** For each of the three bodies recommended in Task 1, produce a native
FreeCAD solid and report, in `src/freecad/geometries/PILOT.md`:

- the method used (RE workbench surface fit, revolve from profile, loft, etc.);
- RMS and max deviation from the source mesh, and the volume error, against the
  acceptance criterion above;
- **wall-clock time spent on that one body** — this is the number that decides
  whether 245 is a week or a year;
- **which parameters, if any, the result actually exposes.** A fitted B-rep with
  no sketch history and no named dimensions is a different representation, not a
  parametric one. Say plainly which it is.

**Done when** `PILOT.md` reports all four items for three bodies, and the solids
are saved under `src/freecad/geometries/pilot/`.

---

### TASK 3 — Rostrocaudal tiling feasibility  *(Claude, not started)*

Whether the T10–T12 block can be tiled rostrally to reach T6–T9, which the
human aim of the R01 requires (DC leads T8/T9, VC lead ventrolateral T6/T7).
Measure, don't assume: vertebral body size trend across the three levels, cord
centreline curvature across the span, root take-off angle by level, and the
seam error a naive translate-and-join would introduce.

---

## Returned results

*Codex: append findings here, newest last. State what you measured and what you
could not. Do not edit sections above this line.*

### 2026-10-05 — Codex: one-vertebra teaching checkpoint

Mohamed redirected the inventory-first task to demonstrate one vertebra before
batch work, requested readable names, allowed assemblies of parts, and asked
for GUI/API teaching pauses. Inventory and the remaining pilots are deferred.
`src/freecad/geometries/pilot/vertebra_01/Vertebra_01_Lesson_01.FCStd` contains
a fitted 18-pole periodic B-spline sketch and editable PartDesign Pad for the
vertebral body, overlaid on `T8-10 - V1-2.STL`. Measured body height:
24.00000057 mm. Pad: one valid closed solid, volume 23368.018732 mm³; a +1 mm
height/recompute test passed and was restored. First checkpoint: approximately
9.5 min interactive elapsed; construction macro: 2.188 s. This is incomplete:
rim fillet, two source channels and posterior reconstruction remain, so no
whole-vertebra RMS/Hausdorff/volume-error acceptance is claimed. See `PILOT.md`
and `pilot/vertebra_01/LESSON_01.md` for the GUI/API lesson and reproducible code.
New assembly identifier `Vertebra_01` preserves source filename/hash metadata;
anatomical level remains unassigned. Source and existing FCStd files unchanged.
Running GUI reports FreeCAD 27.1 development. Strict source connectivity found
eight boundary edges; two boundary-vertex joins in a separate working copy
(maximum move 0.000010790 mm) give a closed reference volume of 35242.215494 mm³.
This is separate from the already-known tissue-interface incompatibility.

### 2026-10-05 — Codex: candidate results and educational cleanup

At Mohamed's request, removed the earlier lesson Markdown, demo macro,
screenshot, result JSON, and saved lesson FCStd; the historical paths in the
preceding entry no longer exist. Open GUI documents were left untouched.
Reconstruction preparation data/scripts and all candidate artifacts remain.
See `src/freecad/geometries/PILOT.md` and
`pilot/vertebra_01_candidate/REPORT.md` relative to that directory.
The candidate assembly has a valid editable Pad → Fillet → Pocket body and
valid static faceted posterior; native fusion failed. Diagnostic mesh-union
source→candidate RMS 0.03066854 mm, reverse RMS 0.05416842 mm, pooled RMS
0.04401578 mm, sampled maximum 1.82762902 mm, volume error -0.123533%.
No certified Hausdorff result. Candidate fails at posterior attachment;
11.18 min measured investigation, not a completed-body estimate. No other
source body reconstructed; no commits or pushes. Next: fit body–pedicle
transition and explicit shared CAD boundary.
