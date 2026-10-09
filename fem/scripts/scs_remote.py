"""Run a mesh or solve stage on a SLURM cluster as if it were a local worker.

RemoteJob has the same signals as async_runner.AsyncProcessRunner, so a macro
keeps its own progress/finished/failed handlers and only the job service
decides where the work runs. One job is a chain of non-blocking steps:

  1. find the remote $HOME;
  2. sync the code and model inputs that provenance hashes, so the cluster
     produces artifacts this workstation accepts;
  3. upload the run's inputs (lead STLs, params.json with LEAD_DIR rewritten,
     and for a solve the mesh artifacts if the mesh was built here);
  4. sbatch fem/hpc/mesh_solve.sbatch with the profile's options;
  5. poll squeue/sacct and stream the job log (PROGRESS lines drive the bar);
  6. on success fetch the published artifacts back into the local run
     directory, rewrite cluster paths, re-check every sha256 and emit
     finished_ok with the worker's own JSON report.

A dropped connection (VPN) never fails the job: the job keeps running on the
cluster and polling resumes when the host is reachable again. Cancel runs
scancel. The service persists `state()` so a job survives a FreeCAD restart
and is fetched when it finishes (resume()).
"""
import json
import os
import re
import shlex
import shutil
import tempfile
import time
from pathlib import Path

from async_runner import QtCore
import config as C
import scs_jobs as J

POLL_MS = 10_000
FETCH = ['mesh.npz', 'mesh.npz.manifest.json', 'preview.vtp', 'preview.vtp.manifest.json',
         'solution.npz', 'solution.npz.manifest.json', 'solution_bg.npz',
         'solution_bg.npz.manifest.json', 'mesh_report.json', 'solve_report.json']
TERMINAL_SLURM = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY',
                  'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE'}


def signature_inputs():
    """Paths (relative to the repository) whose content provenance hashes,
    plus the cluster job scripts."""
    root = Path(C.ROOT)
    paths = [Path(C.TISSUE_MAP_YAML), root / 'fem' / 'requirements.txt']
    paths += [Path(p) for group in C.TISSUE_BODIES.values() for p in group]
    paths += sorted((root / 'fem' / 'scripts').glob('*.py'))
    paths += sorted((root / 'fem' / 'hpc').glob('*'))
    return sorted({str(p.resolve().relative_to(root)) for p in paths if p.is_file()})


class RemoteJob(QtCore.QObject):
    tick = QtCore.Signal(float)
    progress = QtCore.Signal(int)
    stage_message = QtCore.Signal(str)
    finished_ok = QtCore.Signal(dict)
    failed = QtCore.Signal(str)
    aborted = QtCore.Signal(str)
    updated = QtCore.Signal()          # state() changed; the service persists it

    # Overridable for tests (a fake host that runs commands locally).
    ssh_program = 'ssh'
    rsync_program = 'rsync'
    remote_root = J.REMOTE_ROOT

    def __init__(self, profile, task, json_prefix, args, parent=None):
        super().__init__(parent)
        self.profile = J.validate_profile(profile)
        self.host = self.profile['host']
        self.task = task
        self.stage = {'Mesh': 'mesh', 'Solve': 'solve'}[task]
        self.json_prefix = json_prefix
        self.args = [str(a) for a in args]
        self.local_run = os.path.abspath(self.args[0])
        rel = os.path.relpath(self.local_run, C.ROOT)
        if rel.startswith('..'):
            raise ValueError('Run directory must be inside the repository: %s' % self.local_run)
        self.rel_run = rel
        self.background = '--background' in self.args
        self.job_id = None
        self.remote_home = None
        self.slurm_state = None
        self.offset = 0
        self.report_line = None
        self.log_path = os.path.join(self.local_run, 'hpc-%s.log' % self.stage)
        self.script = os.path.join(C.ROOT, 'fem', 'hpc', 'mesh_solve.sbatch')
        self.process = None            # the service reads this for a local pid
        self._last_progress = -1
        self._running = False
        self._cancelling = False
        self._proc = None
        self._started = None
        self._tail = ''
        self._line = ''
        self._ssh_failures = 0
        self.poll_timer = QtCore.QTimer(self)
        self.poll_timer.setSingleShot(True)
        self.poll_timer.timeout.connect(self._poll)
        self.tick_timer = QtCore.QTimer(self)
        self.tick_timer.timeout.connect(lambda: self.tick.emit(time.time() - self._started))

    # -- public ----------------------------------------------------------
    def is_running(self):
        return self._running

    def state(self):
        return dict(host=self.host, job_id=self.job_id, remote_home=self.remote_home,
                    rel_run=self.rel_run, stage=self.stage, background=self.background,
                    slurm_state=self.slurm_state, offset=self.offset, log_path=self.log_path,
                    json_prefix=self.json_prefix)

    def start(self):
        self._begin()
        self._log('Submitting %s to %s (%s)\n' % (self.stage, self.host, self.profile['sbatch_options']))
        self._ssh('echo "$HOME"', self._have_home, 'find the remote home directory')

    def resume(self, saved):
        """Reattach to a submitted job after a restart: poll, then fetch."""
        self.job_id, self.remote_home = saved['job_id'], saved['remote_home']
        self.offset = int(saved.get('offset', 0))
        self.slurm_state = saved.get('slurm_state')
        self._begin()
        self._log('\nResumed watching SLURM job %s after a restart\n' % self.job_id)
        self._poll()

    def abort(self):
        if not self._running or self._cancelling:
            return
        self._cancelling = True
        if self.job_id is None:
            # Not submitted yet: stop the chain; nothing exists on the cluster.
            if self._proc is not None:
                self._proc.kill()
            self._end()
            self.aborted.emit('cancelled before submission')
            return
        self.stage_message.emit('HPC: cancelling SLURM job %s' % self.job_id)
        self._ssh('scancel %s' % self.job_id, lambda out: self._poll(), 'cancel the job')

    # -- chain -----------------------------------------------------------
    def _begin(self):
        self._running = True
        self._started = time.time()
        self.tick_timer.start(1000)
        os.makedirs(self.local_run, exist_ok=True)
        if not os.path.exists(self.log_path):
            Path(self.log_path).write_text('')

    def _end(self):
        for path in getattr(self, '_tmp', []):
            shutil.rmtree(path, ignore_errors=True) if os.path.isdir(path) else os.unlink(path)
        self._tmp = []
        self._running = False
        self.poll_timer.stop()
        self.tick_timer.stop()
        self.updated.emit()

    def _remote(self, *parts):
        return '/'.join([self.remote_home, self.remote_root] + [p for p in parts if p])

    def _have_home(self, out):
        self.remote_home = out.strip().splitlines()[-1]
        self.stage_message.emit('HPC: syncing code and model inputs to %s' % self.host)
        listing = tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False)
        listing.write('\n'.join(signature_inputs()) + '\n')
        listing.close()
        self._tmp = [listing.name]
        self._ssh('mkdir -p %s' % shlex.quote(self._remote()), lambda _:
                  self._rsync(['--checksum', '--files-from=' + listing.name, C.ROOT + '/',
                               '%s:%s/' % (self.host, self._remote())], self._upload,
                              'sync the code mirror'), 'create the code mirror')

    def _upload(self, out):
        self.stage_message.emit('HPC: uploading run inputs')
        staging = tempfile.mkdtemp(prefix='scs-upload-')
        self._tmp.append(staging)
        params = json.loads(Path(self.local_run, 'params.json').read_text())
        if params.get('LEAD_DIR'):
            params['LEAD_DIR'] = self._remote(self.rel_run, 'lead')
        Path(staging, 'params.json').write_text(json.dumps(params, indent=2))
        include = ['lead/***']
        if self.stage == 'solve':
            include += ['mesh.npz', 'mesh.npz.manifest.json', 'preview.vtp',
                        'preview.vtp.manifest.json', 'mesh_report.json']
        filters = ['--include=' + p for p in include] + ['--exclude=*']
        target = '%s:%s/' % (self.host, self._remote(self.rel_run))
        self._ssh('mkdir -p %s' % shlex.quote(self._remote(self.rel_run)), lambda _:
                  self._rsync(filters + [self.local_run + '/', target], lambda _:
                              self._rsync([staging + '/params.json', target], self._submit,
                                          'upload params.json'), 'upload run inputs'),
                  'create the remote run directory')

    def _submit(self, out):
        self.stage_message.emit('HPC: submitting')
        run = self._remote(self.rel_run)
        command = ['sbatch', '--parsable', '-J', 'scs-%s' % self.stage,
                   '-o', run + '/slurm-%j.log'] + J.sbatch_options(self.profile)
        command += [self._remote('fem', 'hpc', 'mesh_solve.sbatch'), run, self.stage,
                    'background' if self.background else 'nobackground']
        self._ssh(' '.join(shlex.quote(c) for c in command), self._submitted, 'submit with sbatch')

    def _submitted(self, out):
        match = re.search(r'^(\d+)', out.strip().splitlines()[-1] if out.strip() else '')
        if not match:
            return self._fail('sbatch did not return a job id:\n' + out)
        self.job_id = match.group(1)
        self.slurm_state = 'PENDING'
        self._log('Submitted SLURM job %s\n' % self.job_id)
        self.stage_message.emit('HPC: SLURM job %s submitted' % self.job_id)
        self.updated.emit()
        self._poll()

    # -- polling ---------------------------------------------------------
    def _poll(self):
        log = self._remote(self.rel_run, 'slurm-%s.log' % self.job_id)
        command = ('squeue -h -j {id} -o "%T|%r" 2>/dev/null; echo "@@SACCT"; '
                   'sacct -n -X -P -j {id} -o State,ExitCode 2>/dev/null | head -1; echo "@@LOG"; '
                   'tail -c +{off} {log} 2>/dev/null; true').format(
                       id=self.job_id, off=self.offset + 1, log=shlex.quote(log))
        self._ssh(command, self._polled, 'poll the job', transient=True)

    def _polled(self, out):
        queue, _, rest = out.partition('@@SACCT\n')
        sacct, _, text = rest.partition('@@LOG\n')
        if text:
            self.offset += len(text.encode('utf-8'))
            self._consume(text)
        queue = queue.strip()
        if queue:
            state, _, reason = queue.partition('|')
            self.slurm_state = state
            if state == 'PENDING':
                self.stage_message.emit('HPC: queued on %s (%s)' % (self.host, reason or 'waiting'))
            elif state == 'RUNNING':
                self.stage_message.emit('HPC: running on %s, job %s' % (self.host, self.job_id))
            self.updated.emit()
            self.poll_timer.start(POLL_MS)
            return
        final = sacct.strip().split('|')
        state = final[0].split()[0] if final and final[0] else ''
        if state not in TERMINAL_SLURM:
            # Left squeue but sacct has not caught up yet: look again shortly.
            self.poll_timer.start(POLL_MS)
            return
        self.slurm_state = state
        self.updated.emit()
        if state == 'CANCELLED':
            self._end()
            self.aborted.emit('SLURM job %s cancelled' % self.job_id)
        elif state != 'COMPLETED':
            self._fail('SLURM job %s ended %s (exit %s)\n%s' % (
                self.job_id, state, final[1] if len(final) > 1 else '?', self._tail[-4000:]))
        elif self.report_line is None:
            self._fail('SLURM job %s completed without a %s line\n%s' % (
                self.job_id, self.json_prefix.strip(), self._tail[-4000:]))
        else:
            self._fetch()

    def _consume(self, text):
        self._log(text)
        self._tail = (self._tail + text)[-64000:]
        self._line += text
        while '\n' in self._line:
            line, self._line = self._line.split('\n', 1)
            if line.startswith(self.json_prefix):
                self.report_line = line[len(self.json_prefix):]
            elif line.startswith(('RESOURCE ', 'Contact ', 'Active:', '  classified ')):
                self.stage_message.emit(line.strip())
            elif line.startswith('PROGRESS '):
                try:
                    pct = max(0, min(100, int(line.split(None, 1)[1])))
                except (IndexError, ValueError):
                    continue
                if pct != self._last_progress:
                    self._last_progress = pct
                    self.progress.emit(pct)

    # -- results ---------------------------------------------------------
    def _fetch(self):
        self.stage_message.emit('HPC: fetching results')
        filters = ['--include=' + name for name in FETCH] + ['--exclude=*']
        self._rsync(filters + ['%s:%s/' % (self.host, self._remote(self.rel_run)),
                               self.local_run + '/'], self._fetched, 'fetch results', transient=True)

    def _fetched(self, out):
        import artifacts
        remote_run = self._remote(self.rel_run)
        for name in ('mesh_report.json', 'solve_report.json'):
            path = Path(self.local_run, name)
            if path.is_file():
                path.write_text(path.read_text().replace(remote_run, self.local_run))
        for manifest in Path(self.local_run).glob('*.manifest.json'):
            target = manifest.with_name(manifest.name[:-len('.manifest.json')])
            try:
                expected = json.loads(manifest.read_text())['sha256']
            except (OSError, ValueError, KeyError):
                continue
            if target.is_file() and artifacts.file_hash(str(target)) != expected:
                return self._fail('%s does not match its manifest after the fetch' % target)
        try:
            report = json.loads(self.report_line.replace(remote_run, self.local_run))
        except ValueError as exc:
            return self._fail('Invalid worker report: %s' % exc)
        report['_stdout_tail'] = self._tail[-4000:]
        report['log_path'] = self.log_path
        report['slurm_job'] = self.job_id
        self._end()
        self.progress.emit(100)
        self.finished_ok.emit(report)

    # -- plumbing --------------------------------------------------------
    def _fail(self, message):
        self._log('\n' + message + '\n')
        self._end()
        self.failed.emit(message + '\nLog: ' + self.log_path)

    def _log(self, text):
        with open(self.log_path, 'a') as f:
            f.write(text)

    def _ssh_options(self):
        control = os.path.join(tempfile.gettempdir(), 'scs-ssh-%C')
        return ['-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', '-o', 'StrictHostKeyChecking=yes',
                '-o', 'ControlMaster=auto', '-o', 'ControlPath=' + control, '-o', 'ControlPersist=120']

    def _ssh(self, command, then, what, transient=False):
        self._start(self.ssh_program, self._ssh_options() + [self.host, 'bash -lc ' + shlex.quote(command)],
                    then, what, transient)

    def _rsync(self, args, then, what, transient=False):
        shell = ' '.join([self.ssh_program] + [shlex.quote(o) for o in self._ssh_options()])
        self._start(self.rsync_program, ['-a', '--partial', '-e', shell] + args, then, what, transient)

    def _start(self, program, args, then, what, transient):
        proc = QtCore.QProcess(self)
        self._proc = proc

        def finished(code, status):
            out = bytes(proc.readAllStandardOutput()).decode('utf-8', 'replace')
            err = bytes(proc.readAllStandardError()).decode('utf-8', 'replace')
            proc.deleteLater()
            if self._proc is proc:
                self._proc = None
            if not self._running:
                return
            if code == 0 and status == QtCore.QProcess.NormalExit:
                self._ssh_failures = 0
                return then(out)
            if transient and self.job_id is not None and code == 255:
                # Connection problem (VPN): the job keeps running remotely.
                self._ssh_failures += 1
                self.stage_message.emit('HPC: %s unreachable (%d tries); job %s keeps running there; retrying'
                                        % (self.host, self._ssh_failures, self.job_id))
                return self.poll_timer.start(POLL_MS) if what == 'poll the job' else \
                    QtCore.QTimer.singleShot(POLL_MS, self._fetch)
            self._fail('Could not %s (exit %s):\n%s' % (what, code, (err or out)[-3000:]))
        proc.finished.connect(finished)
        proc.errorOccurred.connect(lambda error: error == QtCore.QProcess.FailedToStart and finished(127, 0))
        proc.start(program, args)
