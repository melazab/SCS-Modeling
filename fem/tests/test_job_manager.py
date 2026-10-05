"""Profiles, lifecycle persistence, and real local workers behind Job Manager."""
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import scs_jobs as J
from async_runner import QtCore, AsyncProcessRunner
import scs_job_service as S


class JobTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app=QtCore.QCoreApplication.instance() or QtCore.QCoreApplication([])

    def test_settings_snapshot_and_transitions(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=J.Store(tmp);p=store.profiles['local']
            job=store.new_job({'name':'a'},'Mesh',p,['input'],{},None)
            p['mesh_threads']=999
            self.assertNotEqual(job['profile']['mesh_threads'],999)
            store.transition(job,'running');store.transition(job,'completed')
            with self.assertRaises(ValueError):store.transition(job,'running')
            self.assertEqual(store.history()[0]['state'],'completed')

    def test_local_caps_and_remote_probe(self):
        p=J.defaults()['local'];p.update(mesh_threads=100000,solve_cpus=100000)
        validated=J.validate_profile(p)
        self.assertLessEqual(validated['solve_cpus'],J.physical_core_limit())
        self.assertLessEqual(validated['mesh_threads'],J.os.cpu_count())
        p=J.defaults()['case-hpc'];args=J.ssh_probe_args(p)
        self.assertIn('BatchMode=yes',args)
        self.assertIn('StrictHostKeyChecking=yes',args)
        p['host']='case-hpc;touch /tmp/no'
        with self.assertRaises(ValueError):J.ssh_probe_args(p)

    def test_real_worker_retry_and_history(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=J.Store(tmp)
            with patch.object(S.J,'Store',return_value=store):service=S.JobService()
            doc=SimpleNamespace(Name='dorsal',Label='Dorsal',FileName='/tmp/model.FCStd',getObject=lambda _:None)
            script=Path(tmp)/'worker.py'
            script.write_text("print('PROGRESS 35')\nprint('RESULT {\"answer\":42}')\n")
            runner=AsyncProcessRunner(None,sys.executable,str(script),'RESULT ')
            for _ in range(2):
                loop=QtCore.QEventLoop()
                def done(*args):loop.quit()
                runner.finished_ok.connect(done);runner.failed.connect(done)
                service.submit(doc,'Solve',runner,[tmp],{})
                timer=QtCore.QTimer();timer.setSingleShot(True);timer.timeout.connect(loop.quit);timer.start(5000)
                loop.exec();timer.stop()
                runner.finished_ok.disconnect(done);runner.failed.disconnect(done)
                self.assertFalse(runner.is_running())
            self.assertEqual(len(service.records),2)
            self.assertTrue(all(r['state']=='completed' for r in service.records.values()))
            self.assertTrue(all(r['log_path'] for r in service.records.values()))
            self.assertEqual(len(store.history()),2)

    def test_manager_cancels_real_worker(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=J.Store(tmp)
            with patch.object(S.J,'Store',return_value=store):service=S.JobService()
            doc=SimpleNamespace(Name='a',Label='A',FileName='',getObject=lambda _:None)
            script=Path(tmp)/'worker.py';script.write_text('import time\ntime.sleep(30)\n')
            runner=AsyncProcessRunner(None,sys.executable,str(script),'RESULT ')
            loop=QtCore.QEventLoop();runner.aborted.connect(lambda _:loop.quit())
            service.submit(doc,'Mesh',runner,[tmp],{})
            ident=next(iter(service.records))
            QtCore.QTimer.singleShot(100,lambda:service.cancel(ident))
            timer=QtCore.QTimer();timer.setSingleShot(True);timer.timeout.connect(loop.quit);timer.start(5000)
            loop.exec();timer.stop()
            self.assertFalse(runner.is_running())
            self.assertEqual(service.records[ident]['state'],'cancelled')

    def test_remote_never_falls_back_to_local(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=J.Store(tmp)
            with patch.object(S.J,'Store',return_value=store):service=S.JobService()
            doc=SimpleNamespace(getObject=lambda _:SimpleNamespace(Profile='case-hpc'))
            with self.assertRaisesRegex(ValueError,'Remote submission'):
                service.local_settings(doc)

if __name__=='__main__':unittest.main()
