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

## HPC profile

The Case HPC profile is an SSH host alias plus one line of sbatch options, the
same ones you would give `sbatch` or `srun` (for example
`-n 24 --mem=140gb --time=06:00:00`). Your default SLURM account and partition
apply unless you add `-A` or `-p`. Profiles saved with the earlier separate
account/partition/CPU/memory/time fields convert automatically.

**Check cluster with these options** is read-only and submits nothing. It
loads the worker environment (`fem/hpc/scs-env.sh`), reports each tool as OK
or MISSING (SLURM commands, Elmer, MPI, the Python venv), and reads the real
limits: the partition's MaxMemPerCPU, the account's group CPU limit, the CPUs
your group is already running, and the partition's idle CPUs (`sinfo`). The
verdict says how many CPUs SLURM will actually charge. A memory request above
MaxMemPerCPU silently raises CPUs per task (`-n 24 --mem=200gb` became 48),
and a job over the group limit pends forever; `sbatch --test-only` reports
neither. It uses your existing SSH keys and a trusted host key; no passwords
are stored.

Submitting from the panel is not built yet. The service rejects a remote
destination instead of silently running locally.

## Running on Pioneer by hand (until submission is built)

Results are ordinary files on disk at both ends. Nothing is held only in
FreeCAD's memory, and nothing comes back on its own.

1. The worker environment is set up on Pioneer: Elmer in `~/opt/elmerfem-26.2`,
   a Python venv in `~/scs/venv`, and a mirror of the repository in
   `~/scs/SCS-Modeling`. `fem/hpc/scs-env.sh` loads all three.
2. Copy changed `fem/scripts/*.py` files to the mirror with rsync. Provenance
   hashes the code, so the mirror must match the workstation byte for byte, or
   the results will not validate here.
3. Make a run directory under `fem/out/lead_runs/` holding `params.json` and
   `lead/`. Submit it with `sbatch -n <tasks> --mem=<GB> --time=<limit>
   fem/hpc/mesh_solve.sbatch <run_dir> both`.
4. Account `tlv` allows 24 CPUs across the whole group. Batch nodes allow at
   most 6 GB per CPU, so 24 tasks can have at most 144 GB.
5. `fem/hpc/fetch_run.sh <run_name>` copies the published artifacts (mesh,
   preview, solution, manifests, reports, lead STLs) into
   `fem/out/lead_runs/<run_name>`. It rewrites the cluster paths in the reports
   and re-checks every sha256. Elmer's scratch directory (`elmer*/`, about 30 GB
   for a 52 M-tet mesh), `mesh.msh` and the size field stay on the cluster.

The Elmer build is not part of the solution signature, so a field solved on
Pioneer validates on the workstation. Every published field has passed the
independent residual check, and that check is what certifies it. See
`fem/scripts/artifacts.py`.

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
