# FEM execution roadmap

## Current implementation

SCS Potential Visualizer and the batch workflow now use Elmer MPI with
Hypre/BoomerAMG-preconditioned CG as their sole production backend. Contact
basis fields follow the actual selected lead. Local CPU/memory controls,
bounded element preparation and independent residual checks are implemented.
The previous Python solve is retained as a regression reference. See
[ELMER_BACKEND.md](ELMER_BACKEND.md) for validation and installation changes.

## Completed local work

- [x] Add an explicit local execution option in Potential Visualizer, with
  configurable CPU/thread limits and memory budget. Distinguish execution
  location from numerical backend (production backend is Elmer).
- [x] Implement bounded-memory assembly for the Python reference: calculate
  element contributions in batches and accumulate into a global sparse
  system without retaining all element matrices/triplets simultaneously.
  Choose batch size from available/allocated memory and worker count.
  Profile assembly and AMG setup separately; a smaller assembly peak does
  not guarantee the global system and AMG hierarchy fit.
- [x] Validate batched assembly against existing matrices, analytic cases,
  residual/current checks, and measured peak RSS before changing memory gates.
- [x] Implement and validate an Elmer MPI backend using the same classified mesh, tissue
  conductivities, contact source distributions, return/gauge conventions,
  units and result schema. Preserve numerical equivalence; mesh partitioning
  alone is not a replacement for the electrical model or validation.

## Post-poster spatial verification and meshing redesign — deferred

Deferred until after the poster at the user's request. The current poster
results remain preliminary; internal material coverage and mesh convergence
have not been established for the current full dorsal/ventral models.

- [ ] Verify spatial tissue classification in both poster meshes. Plot material
  labels and background elements on transverse/longitudinal sections through
  the cord, CSF, dura, roots and electrodes. Distinguish intended surrounding
  fill from unexpected internal unassigned regions; inspect connected regions
  and quantify their volumes and locations. Connectivity to the exterior alone
  does not establish that a region is anatomically appropriate background.
- [ ] Check tissue coverage, gaps, overlaps and classification precedence against
  the source anatomy using independent containment checks where possible.
  Inspect boundary-crossing tetrahedra and thin-layer resolution. Individually
  watertight surfaces and zero direct epidural–CSF leakage do not establish
  correct material assignment everywhere. Repair or anatomically justify any
  internal gaps instead of automatically assigning background conductivity.
- [ ] Reassess the meshing method before large HPC refinement runs: evaluate
  geometry repair/partitioning and a mesh that follows shared tissue interfaces
  against the current box mesh with centroid-based material assignment. Finer
  elements alone do not repair missing anatomy or incorrect tissue labels.
- [ ] Redesign mesh controls around editable anatomical sizing groups (fine,
  medium, coarse, and lead components), with per-tissue/interface overrides,
  explicit surrounding-domain size and growth controls, and saved comparison
  presets. Keep material definitions separate from mesh sizing groups.
- [ ] Perform mesh-convergence and background-conductivity/domain-size
  sensitivity studies on HPC using predefined potential/E-field metrics in
  matched anatomical regions. Record quality, material volumes, runtime and
  peak memory; select refinement from convergence rather than total tet count.
- [ ] Profile size-field preparation, native Gmsh meshing, classification,
  preview/export, Elmer setup/solve per contact, and result verification/I/O.
  Optimize measured bottlenecks; assess matrix/preconditioner reuse across
  contact basis solves before considering a controller-language rewrite.
- [ ] Benchmark one fixed contact problem at 4, 8 and 14 MPI processes with
  identical solver tolerances and independent verification. Record CPU affinity,
  setup/solve/I/O time and peak memory; do not infer speedup from CPU allocation.
  The completed dorsal and ventral 16-basis workflows took 9,497.2 and 9,073.4 s
  locally. See [measured poster runs](../docs/poster_fast_workflow.md).

## Planned remote work — not implemented

- [ ] Add remote execution over SSH with SLURM. Keep this orchestration out
  of numerical assembly and out of Mesh Generator's geometry responsibilities.
  Configure host/SSH profile, remote work directory, software environment,
  account, partition/QOS, wall time, nodes, MPI tasks, CPUs per task and memory.
  Use existing SSH authentication; do not store credentials in CAD documents.
- [ ] Package immutable mesh/material/lead inputs, solver version and manifests;
  transfer to a job-specific directory. Submit via `sbatch`, launch compute
  steps using the site's MPI/`srun` requirements, record job IDs, expose queued,
  running, failed and completed states, logs, cancellation and reconnection.
  An SSH disconnect or closed FreeCAD window must not lose job tracking.
- [ ] Retrieve results atomically and validate provenance, contact count,
  residuals and completeness before plotting. A changed live lead must not
  silently receive results from a previous submission.
- [ ] Configure remote Elmer/MPI/Hypre software profiles and site-specific
  launch commands. Use the existing distributed backend for multi-node jobs;
  launch the Python preparation/controller once, not one copy per MPI rank.
  Validate CPU affinity, partition counts and process placement through SLURM.

## Resource model

Batched assembly is useful locally and on a remote node. In a distributed
solver, each MPI rank can also batch its owned elements. Batching bounds
temporary memory; it does not itself parallelize the global solve. More CPUs
require an implementation that uses them, and concurrent batches increase
aggregate memory. On one host, processes share its physical RAM. Across nodes,
usable aggregate RAM requires a distributed solver.

References: [Elmer parallel solver documentation](https://www.nic.funet.fi/index/elmer/doc/ElmerSolverManual.pdf)
and [SLURM resource requests](https://slurm.schedmd.com/sbatch.html).

## Mesh tissue labels disagree with STL containment at tissue boundaries

Isolated tets carry a tissue label that contradicts STL ground truth. Confirmed
case: point (57.265, 82.146, 136.307) is inside the white-matter STL and outside
the grey STL, yet its tet is labelled `grey`. Both neighbouring compartments
along the axon are `white`, and a +/-0.05-0.10 mm perturbation in any direction
returns `white` - an isolated boundary sliver, not a region.

Impact today is limited: the label is a provenance check on axon placement, not
an input to the solve, and phi is P1-interpolated from node values the sliver
shares with its white neighbours. So a mislabelled tet perturbs no sampled
number. It does trip `sample_fem.py`'s tissue guard, which is why
`src/neuron/trajectories/vz_lat1_depth1.2.json` carries an explicit, counted
`max_label_exceptions: 1` with the verification recorded in the file.

Impact if left: any per-tissue postprocessing that trusts labels (per-tissue
field statistics, grey-vs-white current budgets, conductivity audits) inherits
the error, and the guard will keep forcing per-trajectory allowances that dilute
its value.

Likely cause: point-in-shell classification at a facet-coincident location,
where white and grey shells share a surface and the containment test resolves a
boundary tet to the wrong side. Candidate fixes: classify by tet centroid with a
consistent tie-break, majority-vote a tet against its neighbours after
classification, or snap coincident white/grey interface facets before labelling.
Add a post-mesh audit that reports isolated tets whose label differs from every
neighbour - that would have caught this without a failed run.
