"""Cluster jobs end to end against a fake host: ssh, rsync and SLURM commands
that run locally, so submission, polling, log streaming, cancellation, a
dropped connection, fetching and resume-after-restart are all exercised."""
import json
import os
from pathlib import Path
import shutil
import stat
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import uuid

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import config as C
import scs_jobs as J
import scs_remote as R
import scs_job_service as S
from async_runner import QtCore, AsyncProcessRunner

FAKE_SSH = r'''#!/bin/bash
# fake ssh: skip options, drop the host, run the command "remotely" (locally).
while [[ $1 == -o ]]; do shift 2; done
shift
[[ -f "$FAKE/offline" ]] && { n=$(cat "$FAKE/offline"); if (( n > 0 )); then echo $((n-1)) > "$FAKE/offline"; exit 255; fi; }
HOME="$FAKE/home" PATH="$FAKE/bin:$PATH" exec bash -c "$*"
'''
FAKE_RSYNC = r'''#!/usr/bin/env python3
import os, subprocess, sys
args, skip = [], False
for a in sys.argv[1:]:
    if skip: skip = False; continue
    if a == '-e': skip = True; continue
    args.append(a[len('fakehost:'):] if a.startswith('fakehost:') else a)
sys.exit(subprocess.call(['rsync'] + args))
'''
FAKE_SBATCH = r'''#!/bin/bash
# fake sbatch --parsable -J name -o pattern <options...> script run stage bg
args=("$@"); n=${#args[@]}
out=""; for ((i=0;i<n;i++)); do [[ ${args[i]} == -o ]] && out=${args[i+1]}; done
script=${args[n-4]}; run=${args[n-3]}; stage=${args[n-2]}; bg=${args[n-1]}
id=$(( $(cat "$FAKE/next_id" 2>/dev/null || echo 100) + 1 )); echo $id > "$FAKE/next_id"
log=${out//%j/$id}
echo "$script $run $stage $bg ${args[*]}" > "$FAKE/jobs/$id.submit"
# While "queued" (PEND seconds) there is no log file yet, as on a real cluster.
( sleep ${PEND:-0}; bash "$FAKE/worker.sh" "$run" "$stage" > "$log" 2>&1; echo $? > "$FAKE/jobs/$id.exit" ) &
echo $! > "$FAKE/jobs/$id.pid"
echo $id
'''
FAKE_SQUEUE = r'''#!/bin/bash
# fake squeue -h -j ID -o fmt
id=$3
[[ -f "$FAKE/jobs/$id.exit" || -f "$FAKE/jobs/$id.cancelled" ]] && exit 0
kill -0 $(cat "$FAKE/jobs/$id.pid") 2>/dev/null || exit 0
ls "$FAKE"/home/scs/SCS-Modeling/fem/out/lead_runs/*/slurm-$id.log >/dev/null 2>&1 && echo "RUNNING|None" || echo "PENDING|Priority"
'''
FAKE_SACCT = r'''#!/bin/bash
# fake sacct -n -X -P -j ID -o fields
id=$5
if [[ -f "$FAKE/jobs/$id.cancelled" ]]; then echo "CANCELLED by 1|0:15"
elif [[ -f "$FAKE/jobs/$id.exit" ]]; then [[ $(cat "$FAKE/jobs/$id.exit") == 0 ]] && echo "COMPLETED|0:0" || echo "FAILED|1:0"
fi
'''
FAKE_SCANCEL = r'''#!/bin/bash
touch "$FAKE/jobs/$1.cancelled"; pkill -P $(cat "$FAKE/jobs/$1.pid") 2>/dev/null; kill $(cat "$FAKE/jobs/$1.pid") 2>/dev/null; exit 0
'''
WORKER = r'''#!/bin/bash
# Stands in for mesh_solve.sbatch: progress, published artifacts, report line.
run=$1; stage=$2
echo "PROGRESS 10"
sleep ${WORKER_SLEEP:-0.5}
head -c 4096 /dev/urandom > "$run/mesh.npz"
sha=$(sha256sum "$run/mesh.npz" | cut -d' ' -f1)
echo "{\"signature\": \"x\", \"sha256\": \"$sha\"}" > "$run/mesh.npz.manifest.json"
echo "{\"mesh_npz\": \"$run/mesh.npz\"}" > "$run/mesh_report.json"
lead=$(python3 -c "import json;print(json.load(open('$run/params.json'))['LEAD_DIR'])")
echo "PROGRESS 80"
echo "MESH_PREVIEW_JSON {\"preview_vtp\": \"$run/preview.vtp\", \"lead_dir\": \"$lead\", \"stage\": \"$stage\"}"
'''


class RemoteJobTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QtCore.QCoreApplication.instance() or QtCore.QCoreApplication([])

    def setUp(self):
        self.fake = Path(tempfile.mkdtemp(prefix='scs-fake-cluster-'))
        for d in ('bin', 'jobs', 'home'):
            (self.fake / d).mkdir()
        for name, text in [('ssh', FAKE_SSH), ('rsync', FAKE_RSYNC), ('sbatch', FAKE_SBATCH),
                           ('squeue', FAKE_SQUEUE), ('sacct', FAKE_SACCT), ('scancel', FAKE_SCANCEL)]:
            path = self.fake / 'bin' / name
            path.write_text(text)
            path.chmod(path.stat().st_mode | stat.S_IEXEC)
        (self.fake / 'worker.sh').write_text(WORKER)
        os.environ['FAKE'] = str(self.fake)
        self.run_dir = Path(C.LEAD_RUNS) / ('_test_remote_' + uuid.uuid4().hex[:8])
        (self.run_dir / 'lead').mkdir(parents=True)
        (self.run_dir / 'lead' / 'SCS Lead Electrode 1.stl').write_text('solid x\nendsolid x\n')
        (self.run_dir / 'params.json').write_text(json.dumps(dict(H_MIN=1, LEAD_DIR=str(self.run_dir / 'lead'))))
        patches = [patch.object(R.RemoteJob, 'ssh_program', str(self.fake / 'bin' / 'ssh')),
                   patch.object(R.RemoteJob, 'rsync_program', str(self.fake / 'bin' / 'rsync')),
                   patch.object(R, 'POLL_MS', 100),
                   # Syncing 87 MB of STLs per test is pointless here; one file proves the path.
                   patch.object(R, 'signature_inputs', lambda: ['fem/hpc/mesh_solve.sbatch'])]
        for p in patches:
            p.start()
            self.addCleanup(p.stop)
        self.profile = dict(name='Fake HPC', backend='slurm', host='fakehost', sbatch_options='-n 2 --mem=4gb')

    def tearDown(self):
        shutil.rmtree(self.run_dir, ignore_errors=True)
        shutil.rmtree(self.fake, ignore_errors=True)

    def wait(self, *signals, timeout=20000):
        loop = QtCore.QEventLoop()
        got = []
        for signal in signals:
            signal.connect(lambda *a, s=signal: (got.append((s, a)), loop.quit()))
        QtCore.QTimer.singleShot(timeout, loop.quit)
        loop.exec()
        return got

    def job(self):
        return R.RemoteJob(self.profile, 'Mesh', 'MESH_PREVIEW_JSON ',
                           [str(self.run_dir), str(self.run_dir / 'params.json')])

    def remote_run(self):
        return self.fake / 'home' / J.REMOTE_ROOT / os.path.relpath(self.run_dir, C.ROOT)

    def test_mesh_job_submits_streams_and_fetches(self):
        job = self.job()
        seen = []
        job.progress.connect(seen.append)
        job.start()
        got = self.wait(job.finished_ok, job.failed)
        self.assertIs(got[0][0], job.finished_ok, got)
        report = got[0][1][0]
        remote = str(self.remote_run())
        # Cluster paths are rewritten to this machine; the worker saw the cluster lead dir.
        self.assertEqual(report['preview_vtp'], str(self.run_dir / 'preview.vtp'))
        self.assertEqual(report['lead_dir'], str(self.run_dir / 'lead'))
        self.assertIn(remote, (self.fake / 'jobs').joinpath('101.submit').read_text())
        self.assertEqual(json.loads((self.remote_run() / 'params.json').read_text())['LEAD_DIR'], remote + '/lead')
        # The artifact came back intact; the local params.json was not overwritten.
        self.assertEqual((self.run_dir / 'mesh.npz').read_bytes(), (self.remote_run() / 'mesh.npz').read_bytes())
        self.assertEqual(json.loads((self.run_dir / 'params.json').read_text())['LEAD_DIR'], str(self.run_dir / 'lead'))
        self.assertIn(str(self.run_dir), (self.run_dir / 'mesh_report.json').read_text())
        self.assertEqual(seen[0], 10)
        self.assertEqual(seen[-1], 100)
        self.assertIn('PROGRESS 80', Path(job.log_path).read_text())
        self.assertFalse(job.is_running())

    def test_queued_job_without_log_yet(self):
        os.environ['PEND'] = '1'
        self.addCleanup(os.environ.pop, 'PEND')
        job = self.job()
        messages = []
        job.stage_message.connect(messages.append)
        job.start()
        got = self.wait(job.finished_ok, job.failed)
        self.assertIs(got[0][0], job.finished_ok, got)
        self.assertTrue(any('queued on fakehost (Priority)' in m for m in messages), messages)

    def test_cancel_runs_scancel(self):
        os.environ['WORKER_SLEEP'] = '30'
        self.addCleanup(os.environ.pop, 'WORKER_SLEEP')
        job = self.job()
        job.start()
        while job.job_id is None:
            self.wait(job.updated, timeout=5000)
        job.abort()
        got = self.wait(job.aborted, job.failed, job.finished_ok)
        self.assertIs(got[0][0], job.aborted, got)
        self.assertTrue((self.fake / 'jobs' / (job.job_id + '.cancelled')).exists())

    def test_dropped_connection_keeps_job_alive(self):
        job = self.job()
        messages = []
        job.stage_message.connect(messages.append)
        job.start()
        while job.job_id is None:
            self.wait(job.updated, timeout=5000)
        (self.fake / 'offline').write_text('3')   # the next three ssh calls fail like a VPN drop
        got = self.wait(job.finished_ok, job.failed)
        self.assertIs(got[0][0], job.finished_ok, got)
        self.assertTrue(any('unreachable' in m for m in messages), messages)

    def test_failed_worker_reports_log_tail(self):
        (self.fake / 'worker.sh').write_text('echo "Traceback: boom"; exit 3\n')
        job = self.job()
        job.start()
        got = self.wait(job.finished_ok, job.failed)
        self.assertIs(got[0][0], job.failed)
        self.assertIn('FAILED', got[0][1][0])
        self.assertIn('boom', got[0][1][0])

    def test_service_drives_panel_runner_and_resumes_after_restart(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = J.Store(tmp)
            store.profiles['hpc'] = self.profile
            with patch.object(S.J, 'Store', return_value=store):
                service = S.JobService()
            doc = SimpleNamespace(Name='d', Label='D', FileName='/tmp/d.FCStd',
                                  getObject=lambda _: SimpleNamespace(Profile='hpc'))
            runner = AsyncProcessRunner(None, sys.executable, 'mesh_preview.py', 'MESH_PREVIEW_JSON ')
            service.submit(doc, 'Mesh', runner, [str(self.run_dir), str(self.run_dir / 'params.json')], {})
            self.assertTrue(runner.is_running())
            got = self.wait(runner.finished_ok, runner.failed)
            self.assertIs(got[0][0], runner.finished_ok, got)
            record = next(iter(service.records.values()))
            self.assertEqual(record['state'], 'completed')
            self.assertEqual(record['remote']['job_id'], '101')

            # A job still running when FreeCAD closed is picked up again by its SLURM id.
            os.environ['WORKER_SLEEP'] = '2'
            self.addCleanup(os.environ.pop, 'WORKER_SLEEP')
            runner2 = AsyncProcessRunner(None, sys.executable, 'mesh_preview.py', 'MESH_PREVIEW_JSON ')
            service.submit(doc, 'Mesh', runner2, [str(self.run_dir), str(self.run_dir / 'params.json')], {})
            while not any((r.get('remote') or {}).get('job_id') == '102' for r in service.records.values()):
                self.wait(service.changed, timeout=5000)
            ident = next(i for i, r in service.records.items() if r['remote']['job_id'] == '102')
            for timer in runner2.findChildren(QtCore.QTimer):
                timer.stop()           # the old session stops watching
            runner2.delegate._running = False
            with patch.object(S.J, 'Store', return_value=J.Store(tmp)):
                restarted = S.JobService()
            self.wait(restarted.changed, timeout=2000)
            job = restarted.runners[ident]
            self.assertTrue(job.is_running())
            self.wait(job.finished_ok, job.failed)
            self.assertEqual(restarted.records[ident]['state'], 'completed')


if __name__ == '__main__':
    unittest.main()
