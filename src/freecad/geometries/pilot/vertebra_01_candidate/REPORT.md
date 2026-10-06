# Vertebra 01 — separate reverse-engineering candidate

2026-10-05. Source: `T8-10 - V1-2.STL`. Anatomical level remains unassigned.
**Interim finding: this hybrid does not meet the reconstruction acceptance criteria.**
It was built and checked in a separate FreeCAD 27.1 GUI on private Xvfb display
`:91`; the user's live GUI and lesson documents were not edited or saved.
No other source body was reconstructed; no files were committed or pushed.

## Native deliverables and method

Open `Vertebra_01_Assembly_Candidate.FCStd` for the two-component assembly, or
`Vertebra_01_Body_Candidate.FCStd` for the editable body alone. `candidate.png`
previews the assembly. Gold is the fitted body; blue-gray is the faceted
posterior retained from the source. The latter has **not** been surface-fitted
or parametrized, and the components overlap.

The editable body continues the existing 18-pole periodic cubic B-spline
sketch → 24.000000570 mm Pad with two new native features: a 1.5 mm Fillet on
both end rims and a through-all Pocket driven by two constrained circles.
The channel radii are 0.499736331 and 0.499731974 mm. Those are neutral source
features; their anatomical purpose has not been assigned.

The posterior is a faceted solid made from the seam-closed source clipped at
local y = 39.5 mm. This geometric cut is **not an anatomical segmentation**.
Unrounded vertices are passed directly to OpenCASCADE; exporting this cut to
binary STL first introduced two zero-area triangles. `prepare_posterior.py`
reproduces the double-precision working data. The source STL is untouched;
its SHA-256 remains `1e80aeff6145e4c938bba5bdfcc0861b766563e294b3a034353a288a63c2d3cf`.

## Validity, parameters and fusion

| Native component | Solids | Faces | Valid / closed | Volume, mm³ |
|---|---:|---:|---|---:|
| Editable body | 1 | 7 | yes / yes | 23220.535781 |
| Faceted posterior | 1 | 25702 | yes / yes | 12236.342621 |

The body exposes Pad Length, Fillet Radius, channel radius and center
constraints, spline poles, and Placement. Height +1 mm, rim radius +0.1 mm,
and one channel radius +0.05 mm each rebuilt to a valid closed single solid;
all were restored. Detailed volumes are in `native_checks.json`.
No anatomical width/depth controls or coupling of posterior anatomy to height
exist. Changing the body does not make the posterior adapt.

Both the normal native CAD fusion and a direct fuzzy fusion with 0.0001 mm
tolerance returned **zero solids and zero faces**. An empty result can still
report `isValid() = True`; the solid-count and positive-volume checks caught
this. The assembly was saved before attempting the failed fusion. There is
no claimed valid native whole-vertebra union or conformal shared interface.

## Whole-exterior diagnostic measurements

Because CAD fusion failed, a separate Manifold mesh Boolean was used strictly
for diagnosis. It is not a replacement native reconstruction. Its in-memory
mesh is watertight and consistently oriented, with 82,876 triangles. The
lossless `diagnostic_mesh_union.npz` preserves this connectivity; the binary
STL export loses exact connectivity at tiny Boolean intersections and must
not be used for topology certification.

The body was tessellated with requested absolute linear deflection 0.005 mm
and angular deflection 0.10 rad. Distances use 100,000 area-uniform samples in
each direction, seed 20261005, and nearest points on target triangles.
These RMS values estimate area averages. The sampled maximum is **not an
exact Hausdorff distance or a certified upper bound**, and no B-rep Hausdorff
certification was performed.

| Diagnostic quantity | Measured result |
|---|---:|
| Source → candidate RMS | 0.03066854 mm |
| Candidate → source RMS | 0.05416842 mm |
| Equal-direction pooled RMS | 0.04401578 mm |
| Bidirectional sampled maximum | 1.82762902 mm |
| Diagnostic union volume | 35198.679556 mm³ |
| Seam-closed source volume | 35242.215494 mm³ |
| Diagnostic relative volume error | -0.123533% |

The sample maximum is far above 0.20 mm; candidate-to-source RMS also exceeds
0.05 mm. Volume agreement alone therefore does not establish accuracy. The
largest discrepancies are at the posterior attachment/cut transition. Most
of the retained posterior is unchanged source geometry, so low global RMS
must not be mistaken for successful posterior reverse engineering.

## Time and next engineering step

Measured wall-clock interval from isolated setup to this report:
**670.5 seconds (11.18 minutes)**.
This includes failed attempts and diagnostic validation, but excludes initial
context reading and the earlier lesson. It is an interim investigation time,
not the time required to finish a vertebra. The latest construction reached
the failed fusion in 69.671 s; independent parameter tests plus the second
fusion attempt took 50.340 s; final distance sampling took
45.811 s. These script runtimes are not substitutes for the
wall-clock investigation time (`timing.json`).

Next, fit the body–pedicle transition from multiple source sections and give
the transition an explicit shared CAD boundary with the body. Separate lofted
or surface-fitted posterior components may be appropriate, but no anatomical
segmentation has yet been approved. Do not scale the constant-profile body
plus arbitrary posterior cut shortcut to the remaining geometries.

Engineering scripts live beside this report; they are not new teaching
macros. `build_candidate.py` and `check_native.py` require an isolated FreeCAD
**GUI**, never `freecadcmd`. Preparation/diagnostics use the temporary Python
environment. No completed Task 2 pilot or patient-specific model is claimed.
