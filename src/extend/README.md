# Extending RADO-SCS rostrally (T6–T9)

`extend_rado.py` adds whole vertebral levels above RADO's top vertebra and
writes a complete, self-contained STL set. `STL_files/` is only read, and
nothing under `fem/scripts/` or `src/ansys/` is modified. Several of those files
are hashed by the solution provenance check; they are only imported.

```bash
fem/.venv/bin/python src/extend/extend_rado.py              # 4 levels: T6-T9
fem/.venv/bin/python src/extend/extend_rado.py --levels 3
fem/.venv/bin/python src/extend/overlap_audit.py fem/out/geometry/rado_up4 --json overlaps.json
```

Output: `fem/out/geometry/rado_up<N>/` (git-ignored, ~180 MB for 4 levels).
It contains RADO's STLs, the new STLs, `manifest.json` and
`extension_check.png`. Ten RADO files are **replaced** in the set, not copied
(listed under `replaced_rado_files`): the 8 continuous bodies, each by one
closed full-length body, and the 2 top sympathetic trunk segments, by joined
versions. The manifest records inputs with SHA-256, every parameter, the
per-level transforms and all check results. A run takes about 6 minutes,
most of it in the checks (`--no-check` skips them).

## Why it works: RADO is one level unit

Measured, and asserted on every run:

- **224 of RADO's 232 per-level bodies are exact rigid copies** of a sibling at
  another level, to about 2e-5 mm (float precision). This covers vertebrae,
  roots, rootlets, DRG, root sleeves, segmental vessels and sympathetic
  ganglia. The exceptions are the 4 discs and one right-side CSF sleeve.
- **But the levels were each placed by hand.** The left DRG steps up 31.3,
  27.6 and 25.8 mm with sideways shifts, while the vertebrae step 31.9 and
  30.8 mm with 4.6° and 5.9° rotations. So the new levels copy **RADO's top
  level arrangement** (top vertebra, the disc above it, and the roots, DRG,
  sleeves, vessels and ganglia at that disc: 63 bodies). They are not an
  average of RADO's levels.
- **Only 8 bodies run continuously:** white, grey, CSF, meninges, epidural,
  two blood layers and the aorta. Each has a constant perpendicular
  cross-section along the cord (Dice ≥ 0.96 over 20 mm).
- **The top 3–8 mm of each continuous body is not the constant section,**
  because the end surfaces cut it. The bodies are therefore cut 10 mm below
  the top. Keep `--join-back` at 10 or more.

## How the levels are stacked

RADO's vertebrae are rigid copies, so the step between two of them is one
rigid motion. The top step, T11 → T10 (disc naming; files `V1-3` → `V1-4`),
is **S: a 5.87° rotation about a left-right axis about 31 cm anterior to the
vertebral bodies**. New level k is S^k applied to the top unit. So:

- **The vertebrae continue RADO's own curve.** The anterior vertebral edge
  of T10–T6 lies on the circular arc that S traces through T11 and T10. RADO
  has no other definition of its curve; it is not a spline. Its lower step,
  T12 → T11, is 4.61°. Four more levels at 5.87° add 23.5° of kyphosis, about
  34° over T6–T12, within normal thoracic kyphosis.
- **Each disc sits on and under its vertebrae exactly as RADO's top disc sits
  on T10.** S carries RADO's T10/T11 disc onto its T9/T10 disc to 0.22 mm, so
  the top disc already fits the T9 vertebra that S makes.
- **The continuous bodies follow the same motion.** Each one's section at the
  cut is carried along the S-orbit (S^t, real t), so every copied vertebra's
  canal, foramina and root entries meet the canal contents exactly as at T10.
  RADO's own cord is not quite S-invariant: at the cut its axis is 1.5° off
  the orbit, so the extension turns by 1.5° there.

An earlier version continued straight, with the same rotation for every
copy. It gave the anterior vertebral edge a sawtooth and left the discs out
of line with the vertebrae above them. Mohamed's review rejected it on
2026-10-12.

## One closed body per continuous tissue

Each continuous body is written as **one closed surface**: the original cut
at a plane, walls swept from the cut points themselves, and a cap at the new
top end. All eight are cut at **one shared plane**, 10 mm below the lowest
top (the cord's), so the nested canal layers keep sections taken at the same
place. Cut each 10 mm below its own furthest point along the cord (up to
1.8 mm apart, see below) they left 2.3× more unassigned canal points.

**The top end.** RADO's canal bodies (grey, white, CSF, dura, epidural and the
two blood layers) end flush in **one plane**, as nested rings that do not
overlap. That plane is tilted about 7° from perpendicular to the cord, so the
wider a body, the further its rim reaches along the cord: up to 1.8 mm more
for the epidural space than for the grey matter. Every body with a cap in that
plane is swept a little past S⁴ of it and cut there, so the extension ends the
way RADO does, flush and four levels higher. The plane is recorded in the
manifest (`rado_top_plane`). The aorta has no cap there and keeps its own end,
12 mm above the cord.

Each top cap is the section's region, holes included: a ring section (an
outer loop and a hole) is closed by one ring of triangles, with each hole
joined to its outline by a slit for the ear clipping.

*Fixed 2026-10-15.* Before then, each body was swept to its own furthest
point. The ends were staggered (white matter 0.25 mm above the grey, CSF
1.0, dura 1.15 and epidural 1.8 mm), which left empty wells up to 0.75 mm deep
over the inner layers. Each loop was also capped on its own: a ring got a full
disc over its outline plus an opposite one over its hole. That is
parity-correct for fem's inside test, but the two coincident discs showed as
dark z-fighting streaks in FreeCAD and left slivers in conform.py's booleans.
Mohamed noticed the top end in FreeCAD.

Cut and walls share vertices, so the result is closed and consistently
oriented. The run verifies this on the file as written (float32): 0 open,
0 non-manifold and 0 misoriented edges for all eight, the same as the RADO
originals. There are no interior caps and no doubled walls, so the same files
serve classification, size fields and display. The earlier "display only"
open surfaces are gone.

## The sympathetic trunk

RADO's trunk segments (`neuro_sympathetic_a`, 29 mm tubes of radius 0.75 mm)
were placed by hand and overlap their neighbours by 0.7–1.1 mm. They are not
S-periodic. S leaves a copied segment 1.1 mm (right) and 1.8 mm (left) away
from the segment below, so a plain copy is disconnected.

Each segment's **upper end is bent** until it overlaps the segment above by
0.8 mm. The bend is a smoothstep from its ganglion to its top end: at most 4.9°
(right) and 10.2° (left), over 15–18 mm. Ganglia and rami stay rigid copies.
RADO's own two top segments get the same bend, so they join the first copy;
they are written as `<stem>__joined.STL`. The topmost copy is left free.

## Overlaps between bodies

RADO's parts were modelled separately and never cut against each other, so
its bodies overlap. `overlap_audit.py` measures this:

| set | different-tissue overlaps | volume |
|---|---|---|
| RADO `STL_files/` | 1,059 pairs | 18.2 cm³ |
| `rado_up4` | 2,079 pairs | 38.8 cm³ |

The copies carry RADO's overlaps unchanged. The FEM labels each tet by tissue
priority (soft tissue before bone), as it always has for RADO. Removing the
overlaps is a requirement for the parametric model. The audit is its
acceptance test.

## Modelling choices

| choice | value | why |
|---|---|---|
| curvature above RADO | RADO's top step, 5.87°/level | continues the curve of RADO's own vertebrae (above) |
| what is copied | RADO's top level unit | levels were placed by hand; the top one is the one adjoining the new levels |
| level names | RADO's disc filenames | Top vertebra = T10, so new levels are T9, T8, T7, T6. The paper's naming is one level higher; see `src/freecad/README.md` "Vertebral levels". New discs are named for the levels they separate (`T_disk_8_9x-1__up1`); other copies keep the source name plus `__up<k>`. Every new name is asserted to map to the same `tissue_map.yaml` class as its source. |
| sizes | unchanged | RADO's own vertebrae are identical copies across T10–T12, and the copies keep that simplification. Real T6–T9 vertebrae are smaller. Size by level belongs to the parametric model. |

## Checks (default run, 4 levels)

Every run checks these on the files it wrote and records them in the
manifest; a failure becomes a warning.

1. **Closure.** Each extended body, read back as written (float32): 0 open,
   0 non-manifold and 0 misoriented edges, the same as the RADO originals.
2. **The cut.** 5 mm below the cut the extended body is the original
   (Dice 1.000). 5 mm above it, the 1.5° turn has moved the section about
   0.1 mm. That gives Dice 0.92–1.00 for thick bodies and 0.71–0.84 for the
   0.5 mm dura and ~1 mm vessels. RADO's own end zone starts there anyway.
3. **Copies.** A section through RADO's top unit, carried by S^k, against
   level k on the same points: per-level tissue agreement ≥ 0.994 and all
   tissue ≥ 0.964, at every level. Sections of the topmost level beyond the
   model's top end are recorded but not gated.
4. **Sympathetic trunk.** 8 of 8 junctions overlap (each lower segment's top
   end lies inside the segment above).
5. **Discs.** Fraction of each disc's vertebra-facing area that touches or
   enters bone: RADO's own discs 0.96–0.99, the new discs 0.96 (below) and
   0.98 (above), identical to RADO's top disc. 95th-percentile gap 0.00 mm
   everywhere.
6. **Gaps.** Points inside the canal outline that belong to no tissue (they
   would become 0.04 S/m background), per section, 0.1 mm grid: extension
   mean 20.2 / max 29, RADO's own canal mean 21.5 / max 86. RADO's own come
   from its mismatched interfaces; the extension adds none beyond them.

`extension_check.png` shows the true midsagittal plane (the arc through all
seven vertebrae) and a curved coronal section along the cord.

## Using the set in the FEM

The `fem/` box mesher labels tets by "inside any body of this class", using
one point-in-shell tester per body. To mesh this set, `fem/` needs a way to
point at a geometry directory instead of `STL_files/`. The 5-entry
`config.TISSUE_STL` (one path per primary compartment, used for size fields,
the cost model and the visualiser surfaces) must then name the extended
bodies (`<stem>__extended_up4.STL`). Each is one closed shell, so one path per
compartment still works. Those files belong to the fem/HPC work and are not
changed here.

Expected cost: the model grows from z ≈ 48–176 mm to about 48–290 mm. At
poster mesh settings, expect ~2.2× the tetrahedra.

A free side result: solving the existing T10 leads on both geometries measures
how much the original model's short length was truncating the field.

## Fused, blended and overlap-free: `conform.py`

```bash
systemd-run --user --scope -q -p MemoryMax=14G -p MemorySwapMax=0 nice -n 19 \
    src/extend/.venv/bin/python src/extend/conform.py --workers 14      # ~50 min
src/extend/.venv/bin/python src/extend/conform_check.py                 # labels vs fem's view of the source
fem/.venv/bin/python src/extend/overlap_audit.py fem/out/geometry/rado_up4_conformal
src/extend/.venv/bin/python src/extend/render_set.py DIR out.png --tissues dura
```

It runs in `src/extend/.venv` (numpy, scipy, manifold3d, scikit-image,
pyvista), not `fem/.venv`, because `fem/requirements.txt` is hashed. Output:
`fem/out/geometry/rado_up4_conformal/`: 89 bodies, a `conform_report.json`, and
per-tissue checkpoints in `stage1/`.

1. **Fuse** each tissue exactly (manifold3d), one cluster of actually
   intersecting bodies at a time. Clustering by bounding box chained all 263
   root bodies into one union and ran out of memory. The sympathetic trunk is
   first rebuilt as one smooth tube. Results: one chain per side, 32 root
   complexes, the dural sac and its sleeves as one body.
2. **Blend** with a 0.25 mm smooth minimum, only where two pieces meet at
   20° or more. The surface is re-made in 4 mm tiles of one global 0.05 mm
   lattice, so neighbouring tiles join without seams. Then simplify to 5 µm.
   **Roots are fused but not blended** (`--no-blend root`). RADO's root bodies
   are layered (a core inside a thin middle layer, rootlets side by side), and
   a blend there adds collars and +0.7% (k 0.1) to +3.4% (k 0.25) root volume
   in the rootlet entry zone.
3. **Partition** by the FEM's priority, piece by piece. Each tissue loses the
   space any higher-priority fused tissue overlapping it claims.

Checked on 2026-10-15, rerun after the top-end fix (see "The top end"):
- 0 open edges in all 89 bodies. 22 have a few pinch edges, where two pieces
  of one body touch after the cut; that does not matter for the inside test.
  The new one is a single 0.0055 mm edge in the white matter, in the top plane
  where its end ring meets the grey matter's.
- `overlap_audit.py`: **0** overlap between different tissues (the source:
  38.8 cm³), 0 within one.
- `conform_check.py` on 400,000 points: labels agree with fem's priority
  labelling of the source set at 99.991% of points spread over the model, and
  99.51% within ~0.15 mm of a surface. No point is in two bodies. The
  differences are the rebuilt trunk and the vessel fillets.
- At the model's top end, the blend's marching cubes bevels the rims of 5 of
  the 7 small epidural vessels (0.55 mm across) by about one 0.05 mm grid step,
  where a blend tile around a vessel junction reaches the cut. Cosmetic.

What the partition shows about RADO, because fem/ already solves it this way:
- Roots lose 29% of their volume to the DRG (dorsal roots run through it) and
  to the cord (rootlets end inside the white matter).
- CSF loses 12% to rootlets and cord-surface vessels.
- Dura loses 22%, mostly the inside of the root sleeves.
- Each vertebra loses a slab of its posterior body, cut off by the epidural
  space, disc and vessels.
- Segmental vessels cut the sympathetic trunk where RADO passes them
  through it.

Interfaces coincide; each side keeps its own triangulation. RADO's gaps
(tissues that never touched) are unchanged.
