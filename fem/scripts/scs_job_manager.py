"""Dedicated FreeCAD Job Manager window. Opening/closing it never changes jobs."""
import json
from pathlib import Path
import time

try:
    from PySide import QtCore, QtWidgets
except ImportError:
    from PySide6 import QtCore, QtWidgets
import scs_jobs as J
from scs_job_service import service


class JobManager(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.service=service()
        self._document=None
        self._editing=False
        self.probe=None
        self.probe_text=None  # shown until the user picks a job row
        outer=QtWidgets.QVBoxLayout(self)
        self.model=QtWidgets.QLabel()
        outer.addWidget(self.model)
        self.destination=QtWidgets.QLabel()
        outer.addWidget(self.destination)
        profiles=QtWidgets.QGroupBox('Execution profiles — settings apply to future jobs')
        form=QtWidgets.QFormLayout(profiles)
        self.profile=QtWidgets.QComboBox()
        for key,p in self.service.store.profiles.items():self.profile.addItem(p['name'],key)
        form.addRow('Edit profile',self.profile)
        self.settings=QtWidgets.QStackedWidget()
        self.local=QtWidgets.QWidget();lf=QtWidgets.QFormLayout(self.local)
        self.mesh_threads=QtWidgets.QSpinBox();self.mesh_threads.setRange(1,J.os.cpu_count() or 1)
        self.cpus=QtWidgets.QSpinBox();self.cpus.setRange(1,J.physical_core_limit())
        self.memory=QtWidgets.QDoubleSpinBox();self.memory.setRange(0,100000);self.memory.setSuffix(' GB');self.memory.setSpecialValueText('Automatic')
        lf.addRow('Mesh maximum threads',self.mesh_threads)
        lf.addRow('Elmer MPI tasks (physical cores)',self.cpus)
        lf.addRow('Solve memory limit',self.memory)
        local_note=QtWidgets.QLabel('Automatic solve memory uses existing worker headroom limits. Meshing retains its own memory estimate. Thread limits do not guarantee full utilization.')
        local_note.setWordWrap(True);lf.addRow(local_note)
        self.settings.addWidget(self.local)
        self.remote=QtWidgets.QWidget();rf=QtWidgets.QFormLayout(self.remote)
        self.host=QtWidgets.QLineEdit();rf.addRow('SSH host',self.host)
        self.sbatch=QtWidgets.QLineEdit();self.sbatch.setPlaceholderText('-n 24 --mem=140gb --time=06:00:00')
        rf.addRow('sbatch options',self.sbatch)
        self.check=QtWidgets.QPushButton('Check cluster with these options')
        self.check.clicked.connect(self.check_connection);rf.addRow(self.check)
        note=QtWidgets.QLabel('Same options you would give sbatch or srun. Your default account and partition apply unless you add -A or -p. '
                              'The check loads the worker environment on the cluster and reads the real SLURM limits; nothing is submitted. '
                              'Submitting from this panel is not built yet: run fem/hpc/mesh_solve.sbatch and bring results back with fem/hpc/fetch_run.sh.')
        note.setWordWrap(True);rf.addRow(note)
        self.settings.addWidget(self.remote);form.addRow(self.settings)
        self.save=QtWidgets.QPushButton('Save profile and use for active model')
        self.save.clicked.connect(self.save_profile);form.addRow(self.save)
        outer.addWidget(profiles)
        self.notice=QtWidgets.QLabel();self.notice.setWordWrap(True);outer.addWidget(self.notice)
        self.table=QtWidgets.QTableWidget(0,5)
        self.table.setHorizontalHeaderLabels(['Model','Task','Destination','State','Progress'])
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self.select_job)
        self.table.cellClicked.connect(self.select_job)
        outer.addWidget(self.table)
        self.cancel=QtWidgets.QPushButton('Cancel selected job')
        self.cancel.clicked.connect(self.cancel_job)
        self.clear=QtWidgets.QPushButton('Clear finished jobs')
        self.clear.setToolTip('Moves completed, failed, cancelled and unverified records to fem/out/jobs/cleared/. Results are not touched.')
        self.clear.clicked.connect(self.clear_finished)
        buttons=QtWidgets.QHBoxLayout();buttons.addWidget(self.cancel);buttons.addWidget(self.clear);outer.addLayout(buttons)
        self.details=QtWidgets.QPlainTextEdit();self.details.setReadOnly(True);outer.addWidget(self.details)
        self.profile.currentIndexChanged.connect(self.load_profile)
        self.load_profile()
        self.timer=QtCore.QTimer(self);self.timer.timeout.connect(lambda: self.refresh() if self.isVisible() else None);self.timer.start(1500)
        self.refresh()

    def load_profile(self):
        p=self.service.store.profiles[self.profile.currentData()]
        remote=p['backend']=='slurm'
        self.settings.setCurrentIndex(int(remote))
        self.save.setText('Save HPC profile' if remote else 'Save profile and use for active model')
        if remote:
            try:p=J.validate_profile(p)  # converts the old account/partition/CPU fields
            except ValueError as exc:self.notice.setText(str(exc))
            self.host.setText(p['host']);self.sbatch.setText(p.get('sbatch_options',''))
        else:
            self.mesh_threads.setValue(p['mesh_threads']);self.cpus.setValue(p['solve_cpus']);self.memory.setValue(p['memory_gb'])

    def edited_profile(self):
        p=dict(self.service.store.profiles[self.profile.currentData()])
        if p['backend']=='local':p.update(mesh_threads=self.mesh_threads.value(),solve_cpus=self.cpus.value(),memory_gb=self.memory.value())
        else:
            p.update(host=self.host.text().strip(),sbatch_options=self.sbatch.text().strip())
        return J.validate_profile(p)

    def save_profile(self):
        import FreeCAD as App
        try:
            p=self.edited_profile();key=self.profile.currentData()
            self.service.store.save_profile(key,p)
            if p['backend']=='local' and App.ActiveDocument:
                self.service.select_profile(App.ActiveDocument,key)
            self.notice.setText('Saved. Running jobs keep their original settings.' if p['backend']=='local' else 'HPC profile saved. Submitting from this panel is not built yet.')
        except (ValueError,OSError) as exc:self.notice.setText(str(exc))
        self.refresh()

    def selected_id(self):
        row=self.table.currentRow()
        return self.table.item(row,0).data(QtCore.Qt.UserRole) if row>=0 and self.table.item(row,0) else None

    def refresh(self):
        import FreeCAD as App
        doc=App.ActiveDocument
        self.model.setText('Active model: '+(doc.Label if doc else 'None'))
        try:self.destination.setText('Next submission: '+(self.service.profile(doc)['name'] if doc else 'Open a model'))
        except ValueError as exc:self.destination.setText(str(exc))
        selected=self.selected_id()
        records=sorted(self.service.records.values(),key=lambda r:r['created'],reverse=True)
        self.table.blockSignals(True);self.table.setRowCount(len(records))
        for row,r in enumerate(records):
            values=[r['document']['label'],r['task'],r['profile']['name'],r['state'],str(r['progress'])+'%']
            for col,value in enumerate(values):
                item=QtWidgets.QTableWidgetItem(value);item.setData(QtCore.Qt.UserRole,r['id']);self.table.setItem(row,col,item)
            if selected==r['id']:self.table.selectRow(row)
        self.table.blockSignals(False)
        self.show_details()

    def select_job(self, *_):
        self.probe_text=None
        self.show_details()

    def set_details(self, text):
        # The 1.5 s refresh must not reset the scroll position or cursor.
        if text==self.details.toPlainText():return
        bar=self.details.verticalScrollBar();at_end=bar.value()==bar.maximum();pos=bar.value()
        self.details.setPlainText(text)
        bar.setValue(bar.maximum() if at_end and pos else pos)

    def show_details(self):
        ident=self.selected_id();runner=self.service.runners.get(ident)
        self.cancel.setEnabled(bool(runner and runner.is_running()))
        self.clear.setEnabled(any(r['state'] in J.TERMINAL for r in self.service.records.values()))
        if self.probe_text is not None:return self.set_details(self.probe_text)
        if not ident:return
        r=self.service.records[ident]
        text=json.dumps(r,indent=2)
        path=r.get('log_path')
        if path and Path(path).is_file():
            with open(path,'rb') as f:
                f.seek(max(0,Path(path).stat().st_size-12000));tail=f.read().decode('utf8','replace')
            text+='\n\nRecent worker log\n'+tail
        self.set_details(text)

    def cancel_job(self):
        ident=self.selected_id()
        if ident:self.service.cancel(ident)
        self.refresh()

    def clear_finished(self):
        self.table.clearSelection();self.details.clear()
        count=self.service.clear_finished()
        self.notice.setText('Cleared %d finished job(s); records moved to %s. Results were not touched.'%(count,self.service.store.root/'cleared'))
        self.refresh()

    def check_connection(self):
        try:self.probe_profile=self.edited_profile();args=J.ssh_probe_args(self.probe_profile)
        except ValueError as exc:self.notice.setText(str(exc));return
        self.probe=QtCore.QProcess(self)
        self.probe.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.probe.finished.connect(self.probe_finished)
        self.probe.errorOccurred.connect(self.probe_error)
        self.check.setEnabled(False);self.notice.setText('Checking the cluster…')
        self.probe.start('ssh',args)
        self.probe_timeout=QtCore.QTimer(self);self.probe_timeout.setSingleShot(True)
        self.probe_timeout.timeout.connect(self.probe.kill);self.probe_timeout.start(20000)
        self.probe_text='Checking the cluster…';self.show_details()

    def probe_error(self,error):
        self.check.setEnabled(True)
        self.notice.setText('SSH check: '+self.probe.errorString())
        self.probe_text='SSH check failed: '+self.probe.errorString();self.show_details()

    def probe_finished(self,code,status):
        self.probe_timeout.stop();self.check.setEnabled(True)
        output=bytes(self.probe.readAllStandardOutput()).decode('utf8','replace')
        verdict,details=J.describe_probe(code,output,self.probe_profile)
        self.notice.setText(verdict)
        self.probe_text=verdict+'\n\n'+details+'\n\nNothing was submitted. Click a job row to return to job details.'
        self.show_details()


def show():
    import FreeCADGui as Gui
    mw=Gui.getMainWindow()
    service().adopt_open_panels()
    dock=mw.findChild(QtWidgets.QDockWidget,'SCS_JobManagerDock')
    if dock is None:
        dock=QtWidgets.QDockWidget('SCS Job Manager',mw)
        dock.setObjectName('SCS_JobManagerDock');dock.setWidget(JobManager())
        mw.addDockWidget(QtCore.Qt.RightDockWidgetArea,dock)
    dock.show();dock.raise_()
    return dock
