# SCS Job Manager — first local implementation

Open `SCS_Job_Manager.FCMacro`. This dedicated window configures execution;
Mesh Generator and Potential Visualizer retain geometry, mesh, and stimulation
controls. Future submissions from the updated macros use the shared profile.
A macro already open before this update keeps its current controls until reopened.
**Do not restart an ongoing solve to obtain the new controls.**

## Available now

- Local profile: mesh maximum threads, physical-core-capped Elmer MPI tasks,
  and solve memory budget (0 = existing automatic headroom behavior).
- Save the local profile and associate it with the active document. Save the
  `.FCStd` to retain its `SCS_Execution.Profile` selection.
- Mesh and solve macros submit through one session service. Profiles are shared;
  editing a profile affects future jobs using it, never jobs already submitted.
- Job table: originating model, task, destination, lifecycle state, progress,
  settings snapshot, command arguments, input metadata where available, worker
  script hash, and recent log. Completed logs are copied to the job archive.
- Explicit cancellation of live jobs, including Elmer child processes through
  the existing subprocess runner.
- **Clear finished jobs** removes completed, failed, cancelled and unverified
  records from the table. Their JSON and archived log move to
  `fem/out/jobs/cleared/`; nothing is deleted, and lead-run results and caches are
  untouched. Running jobs are never cleared.
- Opening Job Manager discovers and monitors already-running macro workers
  without replacing them. These are labeled as existing local jobs.
- Closing/reopening the Job Manager window does not stop jobs. Its service lives
  independently of the dock. The original modeling panel handles result loading.
- Execution references live in `fem/out/jobs/`: `profiles.json`, one job JSON
  per submission, and completed log archives. Retain this directory for history.

## Memory budgets and measured local performance

The profile's explicit memory limit applies to FEM solves. Mesh generation
independently compares its upper peak estimate with **80% of currently available
RAM**. For example, 20 GB available means a 16 GB mesh budget; setting the solve
limit to 22 GB does not increase it. The solve also caps its requested budget at
available-memory headroom. Neither setting reserves physical RAM.

Unused solution caches are released when documents close, without deleting saved
results. The visible mesh panel refreshes available-memory estimates. Do not
change resolution merely to work around a stale estimate.

Both poster models finished locally with 14 MPI processes: dorsal 2 h 38 min,
ventral 2 h 31 min for all 16 bases. Each basis uses those processes sequentially;
14 tasks does not mean 14 contact solutions run concurrently. No core-scaling
benchmark has established the fastest allocation yet. See the
[measured workflow](poster_fast_workflow.md) and [FEM roadmap](../fem/TODO.md).

## State and lifecycle

Local jobs: `starting → running → completed/failed`, or
`running → cancelling → cancelled`. A new submission creates a new record.

The job's profile and metadata are snapshots. Switching documents does not
retarget it. After restarting FreeCAD, unfinished history is marked **unverified**;
we do not assume a reused PID identifies the original process, kill it, or claim
successful recovery. Check its log and cached artifacts. Persistent independent
local-worker recovery and automatic result reattachment after closing a document
remain future work. Keep the originating model open during this first version.

## HPC profile: setup only

The Case HPC profile stores SSH alias, SLURM account, partition, working directory,
CPUs, memory, and wall-time limit. These are profile fields, not model physics.
Blank account/partition/directory fields are allowed during setup.

**Check SSH / SLURM connection** performs a read-only asynchronous probe of host
access and command availability. It submits no jobs, uploads no data, and stores
no credentials. It uses the existing SSH keys/agent and trusted host keys.
VPN or interactive authentication requirements must be satisfied separately.
Missing tools may require loading site-specific environment modules.

Remote submission is intentionally not selectable yet. The service rejects an
unsupported remote destination instead of silently running locally. A successful
SSH probe is not proof that Elmer/Gmsh dependencies or a SLURM allocation work.

## Next remote implementation

1. Discover the user's permitted account/partitions and choose remote storage.
2. Set up and verify the worker environment (Python libraries, Gmsh, Elmer/MPI).
3. Package immutable geometry, metadata, and code; transfer with checksums.
4. Implement SLURM submission, queue polling, reconnect, and explicit cancellation.
5. Run meshing on one high-memory node and Elmer within that allocation initially.
6. Retrieve and validate results, remap remote paths, and associate them only with
   the matching document/model fingerprint.
7. Add lightweight surface/slice/basis retrieval so HPC-sized volume solutions
   do not have to fit in workstation RAM.

Reference inspiration: [ANSYS Workbench HPC Platform Services](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/wb2_help/wb2_help_hpc_platform.html),
particularly separating execution configuration, job monitoring, and result import.
This implementation does not depend on ANSYS services.
