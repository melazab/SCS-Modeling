# Extending RADO-SCS rostrally (T6–T9)

`extend_rado.py` adds whole vertebral levels above RADO's top vertebra and
writes a complete, self-contained STL set. `STL_files/` is only read, and
nothing under `fem/scripts/` or `src/ansys/` is modified. Several of those files
are hashed by the solution provenance check; they are only imported.

```bash
fem/.venv/bin/python src/extend/extend_rado.py                      # 4 levels, straight: T6-T9
fem/.venv/bin/python src/extend/extend_rado.py --levels 3 --kyphosis-deg-per-level 2
```

Output: `fem/out/geometry/rado_up<N>_kyph<deg>/` (git-ignored, ~180 MB for 4
levels). It contains every source STL copied unchanged, the new STLs,
`manifest.json` and `extension_check.png`. The manifest records inputs with
SHA-256, every parameter, the per-level transforms and all check results.
Runtime is about 5 minutes, most of it in the checks (`--no-check` skips them).

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
  cross-section along the cord: Dice ≥ 0.96 over 20 mm and centroid drift
  ≤ 0.03 mm (aorta 0.27 mm). Each is extended by sweeping that section along
  the extended centreline, overlapping the original by 10 mm.
- **The top 3–8 mm of each continuous body is not the constant section,**
  because the end surfaces cut it. The section is therefore taken 10 mm below
  the top. Keep `--join-back` at 10 or more.

Copies and sweeps hang off one centreline frame fitted to the cord: a cubic
y(z), residual 0.017 mm rms. That frame reproduces RADO's own
vertebra-to-vertebra steps to 0.07 mm, which sets the level spacing (30.8 mm
of arc).

## Modelling choices

| choice | default | why |
|---|---|---|
| curvature above RADO | `--kyphosis-deg-per-level 0` (straight) | RADO's cord tilts 17° from z at the bottom and ~1° at the top, so its kyphosis apex is at its *top* vertebra. Repeating RADO's own 5.87° step would bend the spine anteriorly above T10, earlier than real anatomy (apex usually mid-thoracic). |
| where copies are pinned | `--attach roots` | The unit was shaped around a curved canal and is placed on a straight one, so it matches exactly only near the pin. Measured offsets, copy vs RADO, in distance to the cord axis: **pinned at the root complex**, roots/rootlets ≤ 0.16 mm and vertebra ≤ 0.88 mm; **pinned at the fitted vertebra attachment**, rootlets up to 3.5 mm. Bone takes the error because a gap fills with background at 0.04 S/m, the vertebra's own value. |
| level names | RADO's disc filenames | Top vertebra = T10, so new levels are T9, T8, T7, T6. The paper's naming is one level higher; see `src/freecad/README.md` "Vertebral levels". New discs are named for the levels they separate (`T_disk_8_9x-1__up1`); other copies keep the source name plus `__up<k>`. Every new name is asserted to map to the same `tissue_map.yaml` class as its source. |
| sizes | unchanged | RADO's own vertebrae are identical copies across T10–T12; the copies keep that simplification. Real T6–T9 vertebrae are smaller. Size by level belongs to the parametric generator (`cad-dev`). |

## Checks (default run, 4 levels, straight)

1. **Seams.** Swept extension vs original body, mid-overlap: Dice 0.92–0.99
   (dura is the low one, being ~0.3 mm thick).
2. **Copies.** A section through RADO's top unit, carried by each level's
   transform, agrees with the copy at ≥ 0.981 on per-level tissue. Continuous
   bodies are excluded from that gate: they follow the local frame, so carrying
   them with the unit's transform shifts them where RADO curves.
3. **Gaps.** Points inside the canal outline that belong to no tissue (they
   would become background): **extension mean 10.6 / max 24 per section,
   against RADO's own canal mean 34.1 / max 134** (0.1 mm grid). The extension
   adds none. RADO's own come from its mismatched interfaces.

## Using the set in the FEM

The `fem/` box mesher labels tets by "inside any body of this class", using
one point-in-shell tester per body. Overlapping bodies of a class are therefore
a union, and the seams need overlap, not conformity. To mesh this set, `fem/`
needs a way to point at a geometry directory instead of `STL_files/`. The
5-entry `config.TISSUE_STL` (one path per primary compartment, used for size
fields, the cost model and the visualiser surfaces) must also accept the
original body plus its sweep. Those files belong to the fem/HPC work and are
not changed here.

**Open display / size-field surfaces.** Each continuous tissue is two closed
bodies (original + `__sweep_up<N>`), and their caps sit inside each other at
the join. Any code that draws or densifies surfaces would see a fake
interface ring at z ≈ 154–164 and refine around it. So the run also writes
`open_surfaces_display_only/<stem>__open_surface_up<N>.STL`: the original
clipped exactly at the join plane, plus the sweep's walls and top end.
Interior cap area near the joins drops from 1,292 mm² (closed pairs) to 0. The
clip and the sweep's first ring come from the same plane-cut points, so 4 of
the 8 surfaces are actually watertight; the rest have gaps ≤ 0.04 mm. **They
are for visualiser surfaces and size-field clouds only. Never classify them.**
They live in a subdirectory because `config.py` classifies every `.stl` at the
top level of the geometry directory, and the closed bodies stay the
classification set.

Expected cost: the model grows from z ≈ 48–176 mm to about 48–300 mm. At
poster mesh settings, expect ~2.2× the tetrahedra.

A free side result: solving the existing T10 leads on both geometries measures
how much the original model's short length was truncating the field.
