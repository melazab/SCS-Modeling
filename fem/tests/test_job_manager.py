"""Profiles, lifecycle persistence, and real local workers behind Job Manager."""
from pathlib import Path
from types import SimpleNamespace
import subprocess
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

    def test_hpc_profile_is_host_plus_sbatch_options(self):
        old=dict(name='Case HPC',backend='slurm',host='case-hpc',account='',partition='',
                 directory='',cpus=16,memory_gb=64.,time_limit='04:00:00')
        self.assertEqual(J.validate_profile(old),dict(name='Case HPC',backend='slurm',host='case-hpc',
                         sbatch_options='-n 16 --mem=64gb --time=04:00:00'))
        for bad in ['-n 4; rm -rf ~','n 4','--wrap="x y"','$(id)','-n 4 `id`']:
            with self.assertRaises(ValueError):J.validate_profile(dict(old,sbatch_options=bad))

    def test_probe_verdict_reports_memory_cpu_inflation(self):
        out='OK host x\nLIMIT partition=batch account=tlv maxmempercpu=6144 groupcpu=24 groupused=0 cpus=1/2/0/3\n'
        def verdict(options):
            return J.describe_probe(0,out,dict(name='h',backend='slurm',host='h',sbatch_options=options))[0]
        self.assertTrue(verdict('-n 24 --mem=144gb').startswith('Ready: 24 CPUs'))
        # Observed on Pioneer: this became 48 CPUs and pended forever.
        self.assertIn('need 48 CPUs',verdict('-n 24 --mem=200gb'))
        self.assertIn('charge 24 CPUs, not 8',verdict('-n 8 --mem=100gb'))
        self.assertIn('need 40 CPUs',verdict('--ntasks=40 --mem 64G'))
        self.assertIn('missing: ElmerGrid',J.describe_probe(1,'MISSING ElmerGrid\n',
                      dict(name='h',backend='slurm',host='h',sbatch_options='-n 1'))[0])
        self.assertIn('VPN',J.describe_probe(255,'',dict(name='h',backend='slurm',host='h',sbatch_options='-n 1'))[0])

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

    def test_clear_finished_archives_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=J.Store(tmp)
            done=store.new_job({'name':'a'},'Solve',store.profiles['local'],[],{},None)
            store.transition(done,'running');store.transition(done,'completed')
            (Path(tmp)/(done['id']+'.log')).write_text('log')
            live=store.new_job({'name':'b'},'Mesh',store.profiles['local'],[],{},None)
            store.transition(live,'running')
            with patch.object(S.J,'Store',return_value=store):service=S.JobService()
            service.runners[live['id']]=SimpleNamespace(is_running=lambda:True)
            self.assertEqual(service.clear_finished(),1)
            self.assertEqual(list(service.records),[live['id']])
            self.assertEqual([r['id'] for r in store.history()],[live['id']])
            self.assertTrue((Path(tmp)/'cleared'/(done['id']+'.json')).is_file())
            self.assertTrue((Path(tmp)/'cleared'/(done['id']+'.log')).is_file())

    def test_probe_output_survives_refresh(self):
        # Widgets need a QApplication; the other tests already own a QCoreApplication.
        script=PROBE_WIDGET_SCRIPT%str(Path(__file__).resolve().parents[1]/'scripts')
        env=dict(J.os.environ,QT_QPA_PLATFORM='offscreen')
        result=subprocess.run([sys.executable,'-c',script],env=env,capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)


PROBE_WIDGET_SCRIPT='''
import sys,tempfile,types
sys.path.insert(0,%r)
sys.modules['FreeCAD']=types.SimpleNamespace(ActiveDocument=None)
try:
    from PySide import QtWidgets
except ImportError:
    from PySide6 import QtWidgets
app=QtWidgets.QApplication([])
import scs_jobs as J,scs_job_service as S,scs_job_manager as M
tmp=tempfile.mkdtemp();store=J.Store(tmp)
job=store.new_job({'name':'a','label':'A','path':''},'Solve',store.profiles['local'],[],{},None)
store.transition(job,'running');store.transition(job,'completed')
S.J.Store=lambda:store;S._SERVICE=None
M.service=S.service
w=M.JobManager();w.timer.stop()
w.table.selectRow(0)
assert job['id'] in w.details.toPlainText()
w.probe_text='SSH exit code: 0\\nhpc8';w.show_details()
w.refresh();w.refresh()
assert w.details.toPlainText().startswith('SSH exit code'),w.details.toPlainText()[:200]
w.select_job()
assert job['id'] in w.details.toPlainText()
assert w.clear.isEnabled()
w.clear_finished()
assert w.table.rowCount()==0 and not w.clear.isEnabled()
'''

if __name__=='__main__':unittest.main()
