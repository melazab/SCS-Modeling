"""Session service: owned by the application, never by the Job Manager window."""
import copy
import json
import hashlib
import shutil
import os
from pathlib import Path
import time

from async_runner import QtCore
import scs_jobs as J


class JobService(QtCore.QObject):
    changed = QtCore.Signal()

    def __init__(self):
        super().__init__()
        self.store=J.Store()
        self.records={r['id']:r for r in self.store.history()}
        self.runners={}
        self.clients={}
        # Historical local processes cannot safely be adopted from a PID alone.
        # A submitted cluster job can: SLURM still knows it by its job id.
        for r in self.records.values():
            if r['state'] in J.TERMINAL:
                continue
            if r['profile'].get('backend')=='slurm' and (r.get('remote') or {}).get('job_id'):
                QtCore.QTimer.singleShot(0,lambda r=r:self.resume_remote(r))
                continue
            r.update(state='unverified', message='Previous session: check the worker log and cached results. No process has been killed.')
            self.store.write(r)

    def profile_key(self, doc):
        obj=doc.getObject('SCS_Execution')
        return getattr(obj,'Profile','local') if obj else 'local'

    def profile(self, doc):
        key=self.profile_key(doc)
        if key not in self.store.profiles:
            raise ValueError('Execution profile is missing; choose one in SCS Job Manager.')
        return J.validate_profile(self.store.profiles[key])

    def select_profile(self, doc, key):
        if key not in self.store.profiles: raise ValueError('Unknown profile')
        obj=doc.getObject('SCS_Execution')
        if obj is None:
            obj=doc.addObject('App::FeaturePython','SCS_Execution')
            obj.Label='SCS execution settings'
            obj.addProperty('App::PropertyString','Profile','Execution')
            obj.setEditorMode('Profile',1)
        obj.Profile=key
        self.changed.emit()

    def local_settings(self, doc):
        p=self.profile(doc)
        if p['backend']!='local':
            raise ValueError('Remote submission is not implemented; select Local in SCS Job Manager.')
        return p

    def submit(self, doc, task, runner, args, environment):
        profile=self.profile(doc)
        if profile['backend']=='slurm':
            return self.submit_remote(doc,task,runner,args,profile)
        self.track(doc,task,runner,profile,args,environment)
        try:
            runner.start(args,env_overrides=environment)
        except (OSError,ValueError) as exc:
            runner.failed.emit('Could not prepare worker: '+str(exc))
            return
        self.update_running(runner)

    def submit_remote(self, doc, task, runner, args, profile):
        """Run the stage on the profile's cluster. The panel's own runner is
        driven by the remote job, so its handlers and Abort work unchanged."""
        import scs_remote
        try:
            job=scs_remote.RemoteJob(profile,task,runner.json_prefix,args)
        except (ValueError,KeyError) as exc:
            runner.failed.emit('Could not prepare cluster job: '+str(exc))
            return
        runner.attach(job)
        job.setParent(runner)
        self.track(doc,task,runner,profile,args,{})
        ident=next(i for i,r in self.runners.items() if r is runner)
        self._watch_remote(ident,job)
        job.start()
        self.update_running(runner)

    def _watch_remote(self, ident, job):
        def save():
            r=self.records.get(ident)
            if r is None:return
            r['remote']=job.state();r['updated']=time.time()
            self.store.write(r);self.changed.emit()
        job.updated.connect(save)
        job.stage_message.connect(lambda text,i=ident:self.note(i,text))

    def note(self, ident, text):
        r=self.records.get(ident)
        if r is None or r['state'] in J.TERMINAL or not text.startswith('HPC: '):return
        r['message']=text[len('HPC: '):];self.store.write(r);self.changed.emit()

    def resume_remote(self, record):
        """Keep watching a cluster job submitted before a restart; fetch its
        results when it ends. The panel picks them up from the run directory."""
        import scs_remote
        ident=record['id']
        try:
            job=scs_remote.RemoteJob(record['profile'],record['task'],record['remote']['json_prefix'],record['args'],self)
        except (ValueError,KeyError) as exc:
            record.update(state='unverified',message='Cannot resume cluster job: %s'%exc);self.store.write(record);return
        self.runners[ident]=job
        self._watch_remote(ident,job)
        job.progress.connect(lambda value,i=ident:self.progress(i,value))
        job.finished_ok.connect(lambda report,i=ident:self.finish(i,'completed',report))
        job.failed.connect(lambda msg,i=ident:self.finish(i,'failed',msg))
        job.aborted.connect(lambda msg,i=ident:self.finish(i,'cancelled',msg))
        job.resume(record['remote'])
        self.changed.emit()

    def track(self, doc, task, runner, profile=None, args=None, environment=None):
        if runner in self.runners.values(): return
        profile=profile or dict(name='Existing local job',backend='local')
        args=args if args is not None else (runner.process.arguments() if runner.process else [])
        r=self.store.new_job(dict(name=doc.Name,label=doc.Label,path=doc.FileName),task,profile,args,environment,runner.log_path)
        ident=r['id'];self.records[ident]=r
        run_dir=Path(args[0]) if args else None
        if run_dir and run_dir.is_dir():
            r['input_metadata']={}
            for name in ('params.json','mesh.npz.manifest.json','lead/contact_map.json'):
                path=run_dir/name
                if path.is_file():
                    try:r['input_metadata'][name]=json.loads(path.read_text())
                    except (OSError,ValueError):pass
        if Path(runner.script).is_file():
            r['worker_sha256']=hashlib.sha256(Path(runner.script).read_bytes()).hexdigest()
        self.store.write(r)
        self.clients[ident]=runner.parent()  # keep completion callbacks alive if a dock is replaced
        runner.setParent(self)
        self.runners[ident]=runner
        runner.progress.connect(lambda value,i=ident:self.progress(i,value))
        runner.finished_ok.connect(lambda report,i=ident:self.finish(i,'completed',report))
        runner.failed.connect(lambda msg,i=ident:self.finish(i,'failed',msg))
        runner.aborted.connect(lambda msg,i=ident:self.finish(i,'cancelled',msg))
        if runner.is_running(): self.update_running(runner)
        self.changed.emit()

    def owns(self, runner):
        return runner in self.runners.values()

    def update_running(self, runner):
        ident=next((i for i,r in self.runners.items() if r is runner),None)
        if ident is None:return
        record=self.records[ident]
        if record['state'] in J.TERMINAL: return
        record['log_path']=runner.log_path
        record['pid']=runner.process.processId() if runner.process else None
        self.store.transition(record,'running',progress=max(0,runner._last_progress))
        self.changed.emit()

    def progress(self, ident, value):
        r=self.records[ident]
        if r['state'] in J.TERMINAL:return
        r.update(progress=value,updated=time.time())
        self.store.write(r);self.changed.emit()

    def finish(self, ident, state, result):
        r=self.records[ident]
        if r['state'] in J.TERMINAL:return
        # Macro Abort buttons can cancel directly, bypassing manager.cancel().
        if state=='cancelled' and r['state']!='cancelling':self.store.transition(r,'cancelling')
        runner=self.runners[ident]
        if state=='cancelled' and runner.is_running():
            state='failed'
        log_path=runner.log_path
        if log_path and Path(log_path).is_file():
            archived=self.store.root/(ident+'.log')
            shutil.copyfile(log_path,archived)
            log_path=str(archived)
        self.store.transition(r,state,result=result,log_path=log_path)
        self.runners.pop(ident,None)
        self.changed.emit()

    def cancel(self, ident):
        runner=self.runners.get(ident)
        if not runner or not runner.is_running():return
        self.store.transition(self.records[ident],'cancelling')
        self.changed.emit();runner.abort()

    def clear_finished(self):
        finished=[i for i,r in self.records.items() if r['state'] in J.TERMINAL and i not in self.runners]
        for ident in finished:
            self.store.clear(ident);self.records.pop(ident)
        self.changed.emit()
        return len(finished)

    def adopt_open_panels(self):
        # Attach signals only; never replace/restart existing workers or panels.
        import FreeCADGui as Gui
        try:
            from PySide import QtWidgets
        except ImportError:
            from PySide6 import QtWidgets
        for widget in Gui.getMainWindow().findChildren(QtWidgets.QWidget):
            doc=getattr(widget,'doc',None)
            if doc is None:continue
            for attr,task in [('mesh_runner','Mesh'),('solve_runner','Solve')]:
                runner=getattr(widget,attr,None)
                if runner and runner.is_running() and not self.owns(runner):
                    self.track(doc,task,runner)


_SERVICE=None

def service():
    global _SERVICE
    if _SERVICE is None:_SERVICE=JobService()
    return _SERVICE
