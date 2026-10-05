# Local Elmer solve

## Production workflow

Potential Visualizer launches `solve_lead.py`. It validates the selected mesh
and lead, classifies the active tissue domain, exports a native Elmer mesh in
metres, partitions it with ElmerGrid, and runs one unit-current basis solve per
actual contact. Elmer StatCurrentSolver uses Hypre/BoomerAMG-preconditioned CG.
`run_all.sh` invokes the same backend through `run_elmer.py`.

There is one production backend. The previous NumPy/SciPy/PyAMG solver remains
available in source for independent regression comparisons. Its sparse
assembly now also processes bounded element batches instead of constructing
all element matrices and index triplets simultaneously.

The electrical model is preserved: tissue conductivities, finite-conductivity
contact bodies, uniform volumetric +1 A injection for each basis, and a -1 A
uniform-area return over the active domain boundary. Balanced montages cancel
that return. A partitionable point boundary pins one voltage reference.
SIF `Target Nodes` is deliberately avoided: those IDs can refer to different
rank-local nodes after mesh partitioning.

Every result is independently checked by evaluating the P1 matrix action in
batches, without a global Python matrix. Relative residual must be <=1e-6,
pin reaction <=1e-6 A per unit-current solve, all nodes must be present, and
shared MPI nodes must agree. Only completed, verified results receive a
solution manifest. There is no silent fallback to another solver.

## Local controls and diagnostics

- **CPUs:** MPI process count; Python numerical preparation and each rank use
  one numerical-library thread. Sequential preparation may not occupy all CPUs.
- **Memory limit:** decimal GB for the combined worker/child-process RSS.
  Automatic uses 80% of currently available memory. An explicit limit is capped
  at that same available-memory headroom. The worker samples every 250 ms and
  stops itself and solver children if the budget is exceeded. This is a sampled
  safeguard, not an OS-enforced reservation or a guarantee against fast spikes.
  RSS summation conservatively counts shared pages in multiple processes.
- Updated macros use the shared Job Manager profile in `fem/out/jobs/profiles.json`;
  `fem/out/solve_settings.json` supplies legacy defaults. Resource settings affect
  execution, not the mesh cache key. Valid basis fields can be reused after
  changing currents or resource settings.
- **Abort** kills the worker's child processes and its process group, including
  MPI ranks. Logs remain on disk. Live phase/resource messages appear in the panel.
- `solve_lead.py.log` contains phase memory records, independent residuals and
  the final report. Detailed Elmer logs are in `elmer/contactN/solver.log`
  (`elmer_bg/` for the background scenario). Mesh partitioning has its own log.

CLI example (use an existing matching run directory):

```sh
fem/.venv/bin/python fem/scripts/solve_lead.py RUN_DIR \
  --lead-dir RUN_DIR/lead --cpus 4 --memory-gb 16
```

## Full poster runs and performance — September 30, 2026

Both 16-contact background-inclusive runs completed. Dorsal: 17,334,566 tets,
9,497.2 s total, 15.576 GB sampled peak RSS, maximum independent residual
7.103e-9. Ventral: 17,322,431 tets, 9,073.4 s, 15.688 GB peak, residual
1.462e-8. Both used 14 MPI processes per contact. Exact run IDs and reports are
in the [poster workflow](../docs/poster_fast_workflow.md).

The dorsal native Elmer wall times sum to 8,000.93 s over 16 sequential contact
solves; Hypre solution-stage times sum to 7,293.422 s within that total. The
remaining approximately 1,496 s covers work outside those native solver runs.
This establishes where elapsed time went, not a measured 14-core speedup.
Python launches the compiled solver; a C++ controller would not automatically
accelerate the dominant numerical solution stage.

A future 4/8/14-process benchmark should hold mesh, contact, solver settings and
verification fixed, and record affinity, memory and stage times. Evaluate reuse
of the matrix/preconditioner across contact bases: the current implementation
starts a fresh Elmer MPI job for each basis. More cores, finer meshes and
alternative languages are not substitutes for profiling and convergence checks.

These equation residuals establish algebraic solution accuracy only. Internal
material coverage, mesh convergence, anisotropic material behavior and clinical
validation remain unresolved; see [TODO.md](TODO.md).

## Installation on this workstation

The existing MPI-enabled build in `/home/mohamed/Repos/elmerfem/build` was
reconfigured with `WITH_Hypre=TRUE`, using Hypre 3.0.0 from Ubuntu packages
extracted under `/opt/elmerfem/deps`. Its `elmersolver` target was rebuilt and
the resulting core library installed into `/opt/elmerfem/lib/elmersolver`.
The existing MPI executable, ElmerGrid and StatCurrent module remain in use.

Hypre, SuperLU_DIST and their extra Scotch/CombBLAS dependencies are private
to that prefix, rather than installed system-wide. Library search paths were
set with the privately extracted `patchelf` so the installed MPI executable
loads without requiring a shell environment change. The backend also provides
Elmer's data-directory variables to child processes.

Preserved recovery files:

- `/opt/elmerfem/deps/install-before-hypre.tar.gz`: previous `bin`, `lib`, `share`.
- `/opt/elmerfem/deps/CMakeCache.before-hypre.txt`: original build configuration.
- `/opt/elmerfem/deps/packages/`: downloaded dependency packages.

For another installation, enable MPI and Hypre when building Elmer and verify
its runtime dependencies. The current backend expects `/opt/elmerfem/bin`;
remote installation/profile selection is part of the SSH/SLURM TODO.

## Validation

- Unit/integration checks cover batched assembly against a known tetrahedron,
  orientation invariance, matrix-free action, direct-reference versus MPI
  four-contact fields, disconnected domains, memory-budget termination,
  subprocess-tree cancellation, and visualizer contact handling.
- A coarse anatomical mesh (233,591 tetrahedra, including background) passed
  all four Elmer/Hypre basis checks. Worst independent relative residual was
  about 2.5e-9; the fields differed from the accepted Python reference by about
  1.34e-7 of its overall voltage range. Sampled process-tree peak was 0.63 GB
  with two MPI processes. This is a numerical implementation check, not a
  mesh-convergence or physiological validation.
- The native ILU/GCR prototype passed small checks but used about 14 GB and
  converged slowly on the dense mesh. It was stopped and replaced with the
  Hypre configuration; it is not a second selectable backend.
- The user's existing dense mesh completed through the production worker in
  `out/lead_runs/b047c359f91480fea3e651e2`: **1,749,794 active nodes and
  9,894,332 active tetrahedra**, background excluded. Four MPI processes,
  16 GB budget, **9.575 GB sampled peak**, **831.8 seconds** total including
  classification, export, solve, verification and publication. All four
  independent relative residuals were between 1.25e-8 and 1.42e-8; maximum
  absolute pin reaction was 1.42e-10 A. Both mesh and solution manifests pass.
- The requested contact-1 -10 mA / contact-3 +10 mA montage passed nodal
  superposition and zero-net-current checks. Dura, white and grey outputs
  contained no NaNs; the existing surface-sampling fallback used nearby solved
  nodes for 192/1780, 557/5385 and 424/4179 vertices respectively, with maximum
  distances 0.186, 0.204 and 0.193 mm (within the existing 0.5 mm tolerance).
  All three field pipelines were loaded successfully through FreeCAD MCP in
  a temporary document. The original document and its current controls were
  restored without saving CAD changes.

## Remaining limitations

SSH/SLURM submission is not implemented. Local MPI processes still share the
machine's physical RAM. Disconnected active domains are rejected rather than
silently adding grounds or changing material connectivity. A mesh that fits
assembly can still exceed its budget during multigrid setup or iteration.
The previous isotropic-material and mesh-convergence limitations still apply.
