# Handoff — read this first

State of the project as of 2026-09-11. Written for whoever picks this up with
no memory of the conversation that produced it.

## The one-paragraph version

The Ansys current-flow solve has never converged, and after six cluster runs we
know why, and it was never what it looked like. RADO's STL compartments are
**not tessellated compatibly at their shared interfaces** — dura and epidural
share 4 vertices out of ~1780 — so no conformal volume mesh can be built by
stitching them. Workbench worked around that by meshing all 245 bodies
separately and gluing them with 1116 bonded contact pairs, and every solver
failure traces to that: floating bodies, pinball, contact conductance, 878
constraint equations. Meanwhile the independent FreeCAD-side route (`fem/`)
**solved the same problem** by labelling a volume instead of stitching surfaces,
and lands on the published numbers. That route is now the main line.

## What is true, with evidence

- **The interfaces are not conformal.** Round vertices to 1e-4 mm and intersect:
  dura/epidural 4 shared of 1780 and 1756; CSF/dura 4 of 2138 and 1780;
  white/CSF 237 of 5384 and 2138. The STLs themselves are watertight, manifold
  and free of degenerate facets — quality is not the problem, compatibility is.
  This also rules out Workbench "shared topology", which needs coincident
  geometry these bodies do not have.
- **The FreeCAD route has a solved field.** `fem/` — gmsh box mesh graded by
  distance fields, tissue assigned per tetrahedron by point-in-shell
  classification, 308 k nodes / 1.89 M tets, then a P1 FEM. Peak |E| in white
  matter 11.5 kV/m against Khadka's 12; peak surface voltage 1.37 kV against
  1.2. Verified against an analytic sphere (0.61% L-inf) and cross-checked
  against CalculiX on the same mesh (0.0002%). `bash fem/run_all.sh` reproduces.
- **The field is already in the shape NEURON needs.** `fem/out/field_grid.npz`
  holds `phi[i]`, the potential per 1 A into contact i, in V/A, so any
  zero-net-current montage is `sum_i I_i*phi[i]` and bipolar 3-5 is
  `phi[2]-phi[4]`. Eight solves buy the whole configuration space. That was the
  one structural decision `docs/neuron_plan.md` says had to be made early.
- **testA already solves.** `src/ansys/testA_mapdl_cylinder` is a small coaxial
  model that produces a real field (-1.0366 to +1.0380 V for +/-1 mA). Its own
  header records that it works *because* the electrode conductivity was
  moderated to 1e4 S/m and the insulator body dropped, removing the 1e11
  contrast. It also has no contact pairs — nested cylinders `VOVLAP`'d and
  `NUMMRG`'d into one conformal mesh. The physics recipe has been sound all
  along.
- **ECC is not the lever.** Three values four orders apart gave bit-identical
  pivots. **Model size is not the lever either** — stripping a third of the
  elements moved the max pivot from 7.13e16 to 7.73e16, slightly worse.
- **`CEDELE,ALL` is not a fix.** The 878 constraint equations are the electrode
  equipotential definitions. Deleting them lets the model solve by ceasing to be
  the model.

## Mistakes I made — do not repeat these

1. **Three cluster runs lost to the same class of bug: patching a deck by text
   substitution and verifying only that the text landed, never that the result
   was valid APDL in the right processor.** 3797286 OOM (2h35m), 3802282 a bare
   `finish` dropped MAPDL to BEGIN so `solve` was silently ignored (5:48),
   3804186 `CPNGEN` issued in `/SOLU` so every coupled set came out one node
   (35:29). **Rule, now in `src/ansys/preflight/README.md`: nothing goes to the
   full model until it has been rehearsed on testA, where being wrong costs
   minutes.** Mohamed's words: *"You should really plan and think about these
   things before they bite you after the solve is done."*
2. **I wrote a new idiom when a working one was already in the repo.** testA has
   been building equipotential electrodes correctly all along with a single
   `CP,NEXT,VOLT,ALL` over a selected node set. I invented per-node `CPNGEN`
   instead and got 878 sets of one node. Read what is already there first.
3. **I chased the electrode for three runs when the mesh was the problem.** The
   boundary-condition plan was reasonable and is still worth having, but it was
   a workaround for conditioning, not the cause.
4. **I verified colours by counting `ShapeAppearance` blobs.** The count matched
   while 97 of 249 bodies were silently unpainted. Mohamed caught it by looking
   at the render. Reading `DiffuseColor` back off the ViewObject is the check
   that works, and it needs a GUI.
5. **I saved a `.FCStd` from `freecadcmd`**, which silently drops
   `GuiDocument.xml` and every appearance blob. Cost two recoveries. Every
   document-writing script now refuses when `FreeCAD.GuiUp` is false.
6. **A conservation check passed a wrong answer.** In the FEM work, a transposed
   inverse Jacobian left the matrix symmetric with zero row sums, so current
   balance read exactly +/-1.000000 A while the field was wrong by 2.2x and
   refinement did not help. Only the analytic comparison caught it. "Current in
   equals current out" is not a correctness test.

## Plan forward

1. **Make the FreeCAD/gmsh route the main line.** It is the only one that
   addresses the cause, and it is what the published workflow does (Simpleware
   labels a volume; it does not stitch surfaces).
2. **Run Elmer.** Now installed at `/opt/elmerfem` (v26.2, built from source —
   the PPA does not cover this Ubuntu release). `export PATH=/opt/elmerfem/bin:$PATH`
   is all it needs; it finds its own modules. `StatCurrentSolve.so` is present
   and is literally this equation. `fem/out/mesh_tagged.msh` was written with one
   physical volume per tissue so it can go straight into
   `ElmerGrid 14 2 mesh_tagged.msh`. A third independent field, with no analogy
   and no new code.
3. **Fix the C0 problem before NEURON thresholds.** A P1 field is only C0, so
   the activating function (second derivative along an axon) is noisy 10-20 mm
   from the contacts, with visible sign flips, and finer *sampling* makes it
   worse. Refine the cord mesh, go quadratic, or smooth V before differencing.
4. **Add back what was stripped.** The radicular arteries at 0.66 S/m sit right
   beside the electrodes and are currently omitted. Anisotropic white matter too
   (0.1432 transverse / 0.6 longitudinal); it is currently isotropic, which is
   the direction that makes grey matter read ~1.9x high against the paper.
5. **Keep Ansys as a validation target, not the primary path.** If the
   conformal route and Ansys ever agree, that is real mutual validation. Do not
   spend more cluster time patching a mesh topology nobody else uses.
6. **NEURON is at chunk 0 of 8.** `docs/neuron_plan.md` has the breakdown;
   chunks 1-2 (MRG axon, analytic point source) need no FEM at all and can start
   any time.

## Subagents

Two long-running workers have built most of what is in `src/freecad/` and
`fem/`. Both are cold-started per task and know nothing of each other.

- **freecad lead designer** — owns `src/freecad/`. Built `apply_colors.py`,
  `apply_labels.py`, `make_tissue_groups.py`, `build_lead_config.py`,
  `rado_drg_lead.py`, `measure_foramen.py`, `lead_designer.FCMacro`,
  `tissue_visibility.FCMacro`. Knows the FreeCAD traps: the versioned macro
  directory, the headless-save colour loss, the name uniquifier stripping
  trailing digits.
- **freecad fem agent** — owns `fem/` and nothing else. Built the gmsh mesh, the
  P1 solver, the CalculiX cross-check and the field export. It is the one that
  found the non-conformal interfaces.

**Keep them apart.** They have collided once: a worker restored `.FCStd` files
from its own pre-session backups and nearly discarded newer committed work. Tell
each explicitly which directory is its own, and to check
`mcp__freecad__list_documents` before touching anything, because Mohamed usually
has FreeCAD open.

## Working with Mohamed

He is an active collaborator, not a recipient. He caught the unpainted-bodies
bug, the false containment failure, and the DRG lead being wrongly constructed —
all by looking at the model himself. Ask him for `sudo` directly rather than
engineering around a missing package; he has said so explicitly. Use
`-o ConnectTimeout=5` on every `ssh` to `case-hpc` so a dropped CWRU VPN fails
fast instead of hanging.

## Hard-won specifics that will bite again

- Never save a `.FCStd` under `freecadcmd`.
- FreeCAD's macro directory is versioned: `~/.local/share/FreeCAD/v26-3/Macro/`.
  Get it from `FreeCAD.getUserAppDataDir()`.
- Qt renders SVG Tiny: no `clipPath`, and a `--` inside an XML comment is
  illegal XML. Both produce a fully transparent icon with no error.
- The MAPDL deck has CRLF endings, so `$`-anchored `sed` silently never matches.
- An exit code of 0 and a short runtime is not evidence a solve happened. Two
  runs "succeeded" having done nothing. Check for the markers.
- `tissue_map.yaml` conductivities are matched to Khadka on purpose. Do not
  change them to make a solver behave. Mohamed has rejected that once already,
  correctly.
