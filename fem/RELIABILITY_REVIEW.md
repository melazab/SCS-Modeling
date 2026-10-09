# Reliability fixes and validation — 2026-09-15

## September 30 scope update

Both full poster geometries have now completed 16-contact background-inclusive
Elmer runs; measured counts/timings are in the
[poster workflow](../docs/poster_fast_workflow.md). The dated test counts and
smaller-mesh results below are historical evidence, not new tests run today.
Internal spatial classification verification and mesh convergence remain
unestablished for the poster models and are deferred in [TODO.md](TODO.md).
A 0% direct epidural–CSF leak check is not proof of correct material coverage.

## September 15: local Elmer/Hypre backend

Potential Visualizer and `run_all.sh` now use one production solve backend:
Elmer MPI with Hypre/BoomerAMG-preconditioned conjugate gradients. Local CPU
and memory controls, process-tree cancellation, stage memory logging, bounded
element preparation and independent matrix-free residual checks are implemented.
The previous Python solver is retained as a regression reference, with batched
sparse assembly. See [ELMER_BACKEND.md](ELMER_BACKEND.md) for installation,
validation, provenance and limitations; older solver descriptions below are
historical. Remote SSH/SLURM remains unimplemented.

31 automated checks pass, including MPI versus a direct reference, analytical
element matrices, resource guards and child-process cancellation. The three
sphere checks also pass (2.071%, 1.139%, 0.609% maximum interior errors).
On the coarse anatomical validation mesh, all four Hypre fields passed true
residual and pin-reaction checks; maximum relative residual was approximately
2.5e-9, and the maximum voltage difference from the accepted Python reference
was approximately 1.34e-7 of its overall voltage range.

The dense production run also completed successfully: 9,894,332 active tets,
1,749,794 active nodes, four MPI processes, 16 GB budget, 9.575 GB sampled
peak RSS and 831.8 s total. Worst independent relative residual was 1.42e-8.
The accepted result is in `out/lead_runs/b047c359f91480fea3e651e2`; its mesh
and field manifests pass. The user's -10/+10 mA contact-1/contact-3 montage
was checked and plotted in a temporary FreeCAD document. See the backend
document for surface-sampling fallback statistics and remaining limitations.

## Incident

The kernel journal at **10:38:40** reports an OOM kill of Python PID 42307:
`anon-rss:27819292kB` (26.53 GiB). The incomplete run is
`out/lead_runs/274444fa101fe7ce/`, containing a volume mesh but no completed
preview/report. Its 7,630,740 nodes and 47,511,833 tetrahedra explain why the
old global face arrays failed after classification (85% progress).

## Implemented

- Bounded, disk-partitioned face matching for dura checks, tissue surfaces and
  solver boundary extraction; chunked centroid creation; compact preview nodes.
  (2026-10-09: replaced by one in-memory face-adjacency array per stage,
  ~110 bytes/tet, because the bucket files cost ~16 ms per append on Pioneer's
  NFS home and each sort of the 52 M-tet fine mesh took about an hour. See
  `fem/scripts/mesh_faces.py`.)
- Full-anatomy/lead bounding box and shared tissue refinement regions in the
  mesher and estimator. Distance searches are limited to their influence radii;
  structured-grid coordinates are generated in chunks.
- Mesh/preview memory estimates use currently available RAM, with an upper
  bound and extra uncertainty allowance. The worker repeats the check. Solver
  assembly has a separate memory check.
- Mesh Generator owns meshing. Potential Visualizer owns solving and background
  selection; it requires a completed matching mesh and hides stale fields.
- Content-keyed caches, atomic NPZ/manifest publication, completed-preview
  validation, per-run lead snapshots and resumable completed volume meshes.
- Persistent subprocess logs, final-output draining, explicit crash/report
  errors, and silence warnings that do not automatically kill numerical work.
- True residual acceptance at 1e-6, stored residuals/contact IDs/tolerances,
  float64 basis storage and explicit disconnected-mesh rejection. AMG setup is
  seeded and CG can continue in 400-step batches up to 1600 iterations without
  weakening acceptance. Solver-only changes invalidate solutions, not meshes.
- Analytic verification fails if any resolution fails. Background exports use
  `_bg` consistently; basis exports support variable contact counts; tiny-mesh
  outside queries return missing values instead of indexing beyond the mesh.

## Executed checks

| Check | Result |
|---|---|
| 20 targeted tests, including real Qt subprocesses | All pass |
| Analytic sphere, three resolutions | All pass; errors 2.071%, 1.139%, 0.609% |
| Original failed-run mesh: classify, dura check, tissue surfaces | Pass; 46,254,300 active tets, 6,278,730 output triangles |
| Large regression peak resident memory | **3.97 GiB**, 562 s total |
| Coarse full-anatomy integration mesh, live four-contact lead | 40,257 nodes / 233,591 tets; all 11 anatomical classes present |
| Completed volume-mesh cache resume | 0.24 s, without Gmsh |
| Old incomplete cache | Rejected |
| Coarse solve, background removed | Correctly rejected: two disconnected components; no solution published |
| Coarse solve, background included | Four accepted basis fields in 31 s; worst true residual 3.24e-7 |
| Four-contact VTU/grid export | Pass; sample grid shape `(4,3,3,3)` |
| FreeCAD MCP preview loading | 16 tissue/lead meshes, 105,828 facets, including DRG and sympathetic chain |
| Balanced 1 mA c1/c2 montage sampling in FreeCAD | Dura, white and grey outputs; no missing samples |
| FreeCAD dock inspection | Mesh dock has no solve controls; visualizer gates missing meshes |

Validation logs are in `out/large_mesh_check.log` and
`out/lead_runs/full_anatomy_validation/integration_*.log`. Test artifacts are
separate from the historical field of record. The FreeCAD inspection used a
temporary document and restored the original active document; no CAD document
was saved.

## Limits

### Dense preview interaction follow-up

Live run `1426a9c0411db505a88a7881` completed in 419 seconds with
2,637,326 nodes, 16,135,646 raw tets and 9,894,332 classified tets.
Its Gmsh log reports one thread (3D meshing itself took 19.9 seconds).
At diagnosis the core slider set library environment variables, but Gmsh's
`General.NumThreads` remained 1 even with `OMP_NUM_THREADS=20`. This is now
fixed: the worker sets `General.NumThreads` and all three
`Mesh.MaxNumThreads*D` limits from the requested count. Classification/export
do not all parallelize. A separate approximately 491k-tet cube took 0.95 s
with one thread versus 0.53 s with a 20-thread limit; HXT chose up to eight
threads for that workload. This is not a full-anatomy speedup benchmark.

The existing anatomical run was adopted unchanged into
`out/lead_runs/bc7fe5149c7b5a035d1f2254` after verifying geometry, sizing and
artifact hashes. Its report explicitly records compatibility reuse across
the thread-only mesher change, original provenance, and original one-thread
generation. Both original artifacts remain intact; no anatomical remesh ran.

The live document had about 2.1 million visible preview triangles. A mouse
pick over dura took 229 ms, versus 0.15 ms with preview `Selectable=False`.
Panel refresh took 24 ms; an unhighlighted 1000x800 image rendered in 133 ms.
MCP dispatch also deferred calls during interaction; those waits must not
be mistaken for measured rendering time. The reported tens-of-seconds UI
delay has not been fully reproduced in an isolated operation.

Mesh Generator now defaults viewport preview picking off, with a checkbox
to restore it. Tree selection/visibility and the full-resolution mesh remain
available. Both checkbox states were checked through FreeCAD MCP, with
visibility preserved and completed artifact validation still passing.
The setting was applied to the open document without saving the CAD file.
The user confirmed rotation and clicks were noticeably faster afterward.

### Visualizer initialization follow-up

An empty `SCS_Montage.SolutionNpz` previously became `None`, which selected
the CLI's historical eight-contact reference field during document recompute.
That caused the four-versus-eight error simply on opening the panel.
The feature and ViewProvider now remain idle without an explicit solution;
contact-count changes unbind the old solution. Reopening the panel clears
the old binding before attaching refreshed callbacks (attachment itself can
call `updateData`) and hides old field pipelines until a matching result is
selected by Solve / Plot Field.

All 24 terminal tests pass, including real Gmsh option/volume checks and
idle/bound montage tests for differing contact counts. Through FreeCAD MCP,
a temporary document passed fresh-open, 4/8/6/0/4 contact transitions,
stale-solution reopening and idle recomputes with field sampling forbidden.
The user's live four-contact montage is now up-to-date and idle, with the
matching mesh available. No field solve or CAD save was performed.

The integration mesh is deliberately coarse: its dura-leak metric is **1.59%**.
It validates workflow and error handling, not physiological accuracy. The
larger domain, new tissue refinement and revised estimates still need a
convergence study and independent solver validation. The old canal-only
CalculiX/Elmer comparisons do not validate this new anatomical domain.

White matter remains isotropic, root “Middle” material identity remains
unresolved, and P1 field smoothness/NEURON coupling remain future work. No
tissue conductivity was changed. The 1e-6 accepted residual is explicitly
separate from CG's requested 1e-11 target; accepted is not a claim that the
stricter target was reached.
