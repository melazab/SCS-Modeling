# Atlas-defined spinal tracts on the RADO cord

RADO ships one undifferentiated white-matter body, so "dorsal column" and
"spinothalamic tract" have no geometry in the model. This directory maps the
PAM50 white-matter atlas (Spinal Cord Toolbox; De Leener 2018, Lévy 2015) onto
RADO's cord. It is post-processing only: conductivities and the mesh are
untouched, and nothing under `fem/scripts/` is modified (several of those
files are hashed by the solution provenance check; they are only imported).

Outputs go to `fem/out/atlas/` (git-ignored). SCT 7.3 is installed in place at
`~/sct_7.3` without a `PATH` edit; call it as `~/sct_7.3/bin/sct_*`.

## How to run

**Interpreter:** `src/atlas/.venv`, a plain venv with numpy, scipy, nibabel and
matplotlib. `fem/.venv` will not do: it has no nibabel. One-time setup:

```bash
python3 -m venv src/atlas/.venv && src/atlas/.venv/bin/pip install numpy scipy nibabel matplotlib scikit-image manifold3d
```

None of these scripts needs FreeCAD or SCT on `PATH`. Only `pam50_level_match.py`
uses SCT at all, and it just reads PAM50's NIfTI files from
`~/sct_7.3/data/PAM50/template/` (the path is set in the script).
All four read RADO's original `STL_files/` through `fem/scripts/config.py` and
`inside.py`, which are imported, never modified. The three that load geometry set
`SCS_GEOMETRY=rado` before importing config, whose default is the fused T6–T12
set, so they do not see the T6–T9 extension.

Run in this order (each takes seconds):

| # | command | arguments | reads | writes (under `fem/out/atlas/`) |
|---|---|---|---|---|
| 1 | `src/atlas/.venv/bin/python src/atlas/voxelize_cord.py` | `--voxel` mm (default 0.25), `--margin` mm around the CSF sac (1.0), `--out` dir | white, grey, CSF STLs; aorta (left/right check) | `vox<voxel>/rado_{cord,gm,wm}_seg.nii.gz`, `rado_synth_t2.nii.gz`, `voxelize_report.json`; stops if the aorta check finds the frame mirrored |
| 2 | `src/atlas/.venv/bin/python src/atlas/make_disc_labels.py` | `--voxel` (0.25), `--naming discs\|paper` (discs), `--dir` | step 1 output, disc STLs | `vox<voxel>/rado_disc_labels.nii.gz`, `disc_labels.json` |
| 3 | `src/atlas/.venv/bin/python src/atlas/plot_voxel_check.py` | `--voxel` (0.25) | steps 1–2; lead STLs of the two poster runs in `fem/out/lead_runs/` (IDs set in the script) | `vox<voxel>/voxel_check.png` |
| 4 | `src/atlas/.venv/bin/python src/atlas/pam50_level_match.py` | `--z` RADO z of the reference section (110) | step 1 at **0.25 mm** (its centreline); RADO white/grey STLs; PAM50 cord, grey, level maps | `pam50_level_match.json`, `pam50_level_match.png` |
| 5 | `src/atlas/.venv/bin/python src/atlas/register_tracts.py` | `--level` PAM50 vertebral level (T10), `--s` arc length of the RADO section on the extended centreline (70) | RADO white/grey STLs, `FIBRE_CENTRELINE`, PAM50 template + atlas; runs `~/sct_7.3/bin/sct_register_multimodal` and `sct_apply_transfo` | `tracts_<level>/tract_map.{npz,png}`, `tracts.json` |
| 6 | `python src/atlas/build_tract_bodies.py` (Pioneer: `sbatch ~/scs/tract_bodies.sbatch`) | `--atlas`, `--geometry` (the extended conformal set), `--out` | step 5 for T6-T12, extended white/grey and vertebra STLs | `tract_bodies/*.stl`, `bodies.json` (read by `src/freecad/SCS_TractVisualizer.FCMacro`) |

Status: steps 1-4 build and check RADO's cord images and compare RADO's
cross-section with PAM50's levels; step 5, run once per level T6-T12,
registers the atlas onto RADO's section, and step 6 turns the seven maps into
3D tract bodies (below).

## Tract registration (step 5)

Each PAM50 vertebral level T6-T12 is registered separately onto RADO's
constant cross-section (Mohamed, 2026-10-21: "register all levels from T6-T12
matching RADO"), so each extended-model vertebra gets its own level's tract
layout on the same cord. RADO's section is always sampled inside RADO's own
STLs (`--s 70` for T6-T9, whose heights lie in the copied levels). Inputs
are synthetic images, white 1 / grey 2, blurred 0.25 mm; SCT runs affine then
SyN (MeanSquares); the sweep that chose these is in the script's PARAM comment.

| PAM50 level | T6 | T7 | T8 | T9 | T10 | T11 | T12 |
|---|---|---|---|---|---|---|---|
| grey Dice before | 0.46 | 0.47 | 0.39 | 0.45 | 0.44 | 0.62 | 0.58 |
| grey Dice after | 0.79 | 0.78 | 0.75 | 0.78 | 0.77 | 0.91 | 0.95 |
| gracilis, mm² (both sides) | 4.7 | 4.8 | 4.9 | 5.8 | 6.4 | 6.4 | 5.9 |
| cuneatus | 3.4 | 3.4 | 3.3 | 2.2 | 1.7 | 1.9 | 2.1 |
| lateral corticospinal | 3.0 | 3.1 | 2.4 | 2.6 | 3.0 | 3.5 | 4.4 |
| spinal lemniscus | 2.8 | 2.6 | 2.6 | 2.8 | 2.6 | 2.1 | 2.1 |

RADO has no image contrast, so only the white/grey outlines are matched: the
tract layout inside the white matter is the atlas prior, renormalised to sum
to 1 in RADO's white matter. RADO's "H" puts most of its grey in the ventral
horn and intermediate zone labels. Cuneatus below about T6 is an atlas
artefact (see caveats).

`build_tract_bodies.py` (step 6, on Pioneer: `sbatch ~/scs/tract_bodies.sbatch`,
~30 s, in `~/scs/atlas-venv` with scikit-image and manifold3d) turns the seven
maps into one body per tract and side, `fem/out/atlas/tract_bodies/*.stl` plus
`bodies.json`; its docstring has the method. What it guarantees, and the
first-build defect each one fixed:

- the bodies fill exactly the model's cord: masked by the extended white + grey
  STLs (the first build filled the whole centreline and ran ~10 mm past the
  caudal end);
- neighbours share their boundaries, no gaps or overlaps (each tract smoothed
  alone came out 5 % short); total mesh volume 7063 mm³ against 7143 mm³ of
  cord, the 1 % being the rounding of the cord's outer edge;
- flat ends flush with the cord's own oblique end planes, trimmed with
  manifold3d (the 0.5 mm grid made them terraces);
- no corner-only contacts (islands relabelled), the cause of the first build's
  non-manifold edges;
- each vertebra's map sits at that vertebra's centroid (`src/freecad/scs_levels.py`),
  not at its label.

The lateral reticulospinal tracts (PAM50 labels 10, 11) win no voxel of RADO's cord
at any level and have no body; every other tract, both sides, has one.

`src/freecad/SCS_TractVisualizer.FCMacro` shows them in FreeCAD: one toggle
puts the bodies (or one slice per vertebra, from that vertebra's own map) in a
"Tract bodies" group, with a tract tree grouped by funiculus.
`SCS_CrossSectionViewer.FCMacro` cuts them, and the anatomy, square to the
cord at any level.

## Frame — checked, do not re-derive

RADO STL mm are **LPS**: +x patient **left**, +y **posterior**, +z
**superior**. +x = left is established two independent ways: the descending
aorta (`Thoracic_aorta-1.STL`) is 10.8 mm toward +x of the cord, and a
right-handed frame with +y posterior and +z superior forces it. NIfTI world is
RAS = (-x, -y, z), so the voxel-to-RAS affine is a 180° rotation about z, not a
reflection. On the poster runs, lead 1 (C1–C8) is on the patient's right and
lead 2 (C9–C16) on the left, each ~2 mm from the midline.

`voxelize_cord.py` reproduces the STL volumes to 0.5 % (grey) and 0.1 % (white)
and asserts the aorta test on every run.

## Vertebral levels — names are not settled

Disc labels follow SCT's convention (value = vertebra below the disc,
C1 = 1), using RADO's disc filenames, which make the vertebrae T10–T12 (labels
17–20, ~29 mm apart). The paper says T9–T11 and the STL filenames say "T8-10";
see `src/freecad/README.md` "Vertebral levels". `--naming paper` shifts by one.

## RADO's cord is one cross-section, swept

Perpendicular sections at z = 70, 80, 110, 137 and 150 mm agree to Dice
≥ 0.987: 8.2 × 4.8 mm, 31.3 mm² cord, 9.4 mm² grey (30 %). The cord tilts
1–17° from z along its length, so z-slices are oblique; use perpendicular
sections for anything shape-related.

Against PAM50 (`pam50_level_match.py`, cords normalised to their own extents):

| PAM50 | grey fraction | grey Dice vs RADO | spinal lemniscus, one side |
|---|---:|---:|---:|
| vertebral T6–T10 (thoracic segments) | 0.17–0.21 | 0.42–0.48 | 2.2–2.5 mm² |
| vertebral T11 | 0.25 | 0.67 | 1.4 mm² |
| vertebral T12 (lumbar enlargement) | 0.37 | 0.72 | 0.8 mm² |
| RADO | 0.30 | — | (no tracts) |

RADO's grey matter is a schematic broad "H" on a flatter cord (LR/AP 1.7 vs
PAM50's 1.3–1.45); by proportion it sits nearest the lumbar enlargement.

**Consequence for registration.** Standard vertebral-level registration would
give RADO's top vertebra (ventral lead, z ≈ 137) a thoracic tract layout and its
bottom one (dorsal lead, z ≈ 80) a lumbar layout, with a spinal lemniscus about
a third the size, on a cord that is geometrically identical at both. Per-tract
dorsal-versus-ventral differences would then partly be atlas-level artefacts.
The alternative is to register a single PAM50 level onto RADO's constant
cross-section and sweep it along the centreline. Which to use is open.

## Atlas caveats

- Lévy 2015 propagated one drawn cross-section along the cord; fasciculus
  cuneatus is present in PAM50 down to T12, where it does not exist
  anatomically. Treat gracilis + cuneatus as "dorsal columns".
- The spinal lemniscus is spinothalamic **and** spinoreticular.
- RADO has no image contrast: registration fits only white/grey outlines, so the
  internal tract layout is entirely the atlas prior. Report grey-matter Dice
  next to every per-tract number.
