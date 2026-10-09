# SCS Job Manager

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

Jobs: `starting → running → completed/failed`, or
`running → cancelling → cancelled`. A new submission creates a new record.
A cluster job is `running` from submission on; the table adds its SLURM state.

The job's profile and metadata are snapshots. Switching documents does not
retarget it. A cluster job is resumed after a restart by its SLURM job id (see
below). After restarting FreeCAD, an unfinished local job is marked **unverified**;
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

**Save profile and use for active model** makes a model's Generate Mesh and
Solve / Plot Field run on the cluster.

## Running a stage on the cluster

With the HPC profile in use, Generate Mesh and Solve / Plot Field behave as
they do locally: the same progress bar, status line and Abort button. The
status line adds the cluster state (`queued on case-hpc (Priority)`,
`running on case-hpc, job 3996452`). `fem/scripts/scs_remote.py` does the work:

1. **Sync the code mirror.** Every file that provenance hashes (the pipeline
   scripts, the tissue map, the anatomy STLs, `requirements.txt`) and
   `fem/hpc/` is copied to `~/scs/SCS-Modeling` with `rsync --checksum`. A
   result made on the cluster then validates here; uncommitted changes go too.
2. **Upload the run inputs.** The lead STLs and `params.json` go up, with
   `LEAD_DIR` pointed at the cluster copy. For a solve, the mesh artifacts go up
   too if the mesh was built on the workstation.
3. **Submit** `fem/hpc/mesh_solve.sbatch` with the profile's sbatch options.
   CPUs and the solve memory budget come from the allocation.
4. **Poll** every 10 s over one reused SSH connection, streaming the job log
   into `<run>/hpc-mesh.log` or `hpc-solve.log`. A dropped connection (VPN) is
   reported and retried; it never fails the job, which keeps running remotely.
5. **Fetch** on completion: mesh, preview, solution, manifests and reports come
   back into the run directory. Cluster paths in the reports are rewritten, and
   every sha256 is checked against its manifest before the panel loads the
   result. Elmer's scratch (`elmer*/`, about 30 GB at 52 M tets), `mesh.msh`
   and the size field stay on the cluster.

Abort or Cancel runs `scancel`. The job record keeps the SLURM job id. If
FreeCAD closes, the job keeps running, and on restart the job service picks it
up again and fetches the results when it ends. Generate Mesh or Solve / Plot
Field then finds them as a completed cache entry.

The Mesh Generator's memory estimate is judged against the profile's `--mem`
allocation instead of this machine's RAM. The Elmer build is not part of the
solution signature: every published field passes the independent residual
check, and that check is what certifies it (`fem/scripts/artifacts.py`).

## Cluster facts (Pioneer)

- The worker environment: Elmer in `~/opt/elmerfem-26.2` (the workstation's
  revision), a Python venv in `~/scs/venv`, and the repository mirror in
  `~/scs/SCS-Modeling`. `fem/hpc/scs-env.sh` loads all three.
- Account `tlv` allows 24 CPUs across the whole group. Batch nodes allow at
  most 6 GB per CPU, so 24 tasks can have at most 144 GB.
- Nodes are shared. Elmer's time per contact on the same 52 M-tet mesh has
  ranged from 265 s to 1000 s depending on what else ran on the node.
- The home directory is NFS: one small write costs about 16 ms against 0.07 ms
  on local disk, so nothing here writes many small scratch files there.

`fem/hpc/fetch_run.sh <run_name>` still copies a run back by hand, for jobs
submitted outside the panel.

## Not built yet

- Lightweight surface/slice/basis retrieval, so HPC-sized volume solutions
  do not have to fit in workstation RAM.
- Choosing which cached cluster result belongs to an already-open document
  without pressing Generate Mesh or Solve / Plot Field again.

Reference inspiration: [ANSYS Workbench HPC Platform Services](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/wb2_help/wb2_help_hpc_platform.html),
particularly separating execution configuration, job monitoring, and result import.
This implementation does not depend on ANSYS services.
