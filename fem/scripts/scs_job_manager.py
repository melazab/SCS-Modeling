"""SCS Job Manager: where jobs run, and what they are doing. Opening, closing or
rebuilding this window never changes a job; the jobs belong to
scs_job_service, which this only reads and asks (cancel, clear).

Layout (2026-10-23, after the Tract visualizer's rules: one control per
state, one status per job, explanations in tooltips):
  Run on          the active model's execution profile, one switch
  <profile> box   its settings (host and sbatch options, or local threads),
                  Check and Save; Save lights up only when something changed
  Jobs            one Status column per job: the stage and its time while
                  running (scs_stages), "Queued (Priority)", "Done · 1 h 02 min"
  Selected job    a short summary and the worker log's tail; the raw record
                  is one click away (Copy record)
"""
import json
import time
from pathlib import Path

try:
    from PySide import QtCore, QtGui, QtWidgets
except ImportError:
    from PySide6 import QtCore, QtGui, QtWidgets
import scs_jobs as J
import scs_stages
import scs_monitor
from scs_job_service import service

# Read-only cluster check (scs_jobs.ssh_probe_args). See check_connection().
PROBE_TIMEOUT_MS = 90000
STATE_COLOUR = {"running": "#a15c00", "queued": "#a15c00", "cancelling": "#a15c00", "completed": "#1b5e20",
                "failed": "#b00020", "cancelled": "#666666", "unverified": "#666666"}
LOCAL_TIPS = {"mesh_threads": "gmsh threads for meshing on this machine",
              "cpus": "Elmer MPI tasks, physical cores",
              "memory": "Solve memory limit; Automatic uses the worker's headroom limits. Meshing keeps its own "
                        "memory estimate. Thread limits do not guarantee full utilisation."}
REMOTE_TIP = ("Same options you would give sbatch or srun. Your default account and partition apply unless you "
              "add -A or -p (e.g. -p smp for more than 6 GB per CPU). With this profile, Generate Mesh and Solve "
              "sync the code to the cluster, submit there and fetch the results back; jobs keep running if "
              "FreeCAD closes and are picked up again on restart.")


def hint(text=""):
    lab = QtWidgets.QLabel(text)
    lab.setWordWrap(True)
    lab.setStyleSheet("color: gray;")
    return lab


class JobManager(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        import sys, os
        sys.path.insert(0, os.path.join(Path(__file__).resolve().parents[2], "src", "freecad"))
        import scs_ui
        self.ui = scs_ui
        self.service = service()
        self.probe = None
        self.probe_text = None
        self.clock = scs_stages.StageClock()
        self.monitor = scs_monitor.monitor()
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(6)

        # -- Run on: the active model's profile ------------------------------
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Run on"))
        store = self.service.store
        self.keys = sorted(store.profiles, key=lambda k: (store.profiles[k]["backend"] != "local", k))
        self.run_on = scs_ui.Segmented([store.profiles[k]["name"].replace(" workstation", "") for k in self.keys])
        self.run_on.setToolTip("Where the active model's next mesh and solve run. Running jobs keep their own settings.")
        self.run_on.changed.connect(self.use_profile)
        row.addWidget(self.run_on, 1)
        outer.addLayout(row)

        # -- the selected profile's settings ----------------------------------
        self.box = QtWidgets.QGroupBox()
        box = QtWidgets.QVBoxLayout(self.box)
        local = QtWidgets.QWidget()
        self.local_page = local
        lf = QtWidgets.QFormLayout(local)
        lf.setContentsMargins(0, 0, 0, 0)
        self.mesh_threads = QtWidgets.QSpinBox()
        self.mesh_threads.setRange(1, J.os.cpu_count() or 1)
        self.cpus = QtWidgets.QSpinBox()
        self.cpus.setRange(1, J.physical_core_limit())
        self.memory = QtWidgets.QDoubleSpinBox()
        self.memory.setRange(0, 100000)
        self.memory.setSuffix(" GB")
        self.memory.setSpecialValueText("Automatic")
        for w, cap, tip in ((self.mesh_threads, "Mesh threads", "mesh_threads"), (self.cpus, "Solve cores", "cpus"),
                            (self.memory, "Solve memory", "memory")):
            w.setToolTip(LOCAL_TIPS[tip])
            lf.addRow(cap, w)
            w.valueChanged.connect(self.mark_edited)
        box.addWidget(local)
        remote = QtWidgets.QWidget()
        self.remote_page = remote
        rf = QtWidgets.QFormLayout(remote)
        rf.setContentsMargins(0, 0, 0, 0)
        self.host = QtWidgets.QLineEdit()
        self.host.setToolTip("SSH host alias, as in ~/.ssh/config")
        self.sbatch = QtWidgets.QLineEdit()
        self.sbatch.setPlaceholderText("-p smp -n 24 --mem=600gb --time=12:00:00")
        self.sbatch.setToolTip(REMOTE_TIP)
        rf.addRow("Host", self.host)
        rf.addRow("sbatch", self.sbatch)
        for w in (self.host, self.sbatch):
            w.textEdited.connect(self.mark_edited)
        box.addWidget(remote)
        buttons = QtWidgets.QHBoxLayout()
        self.check = QtWidgets.QPushButton("Check cluster")
        self.check.setToolTip("Loads the worker environment on the cluster and reads the real SLURM limits for these "
                              "options (partition memory per CPU, group CPU cap). Nothing is submitted.")
        self.check.clicked.connect(self.check_connection)
        self.save = QtWidgets.QPushButton("Save")
        self.save.setToolTip("Save these settings; they apply to the next job. Running jobs keep their own.")
        self.save.clicked.connect(self.save_profile)
        buttons.addWidget(self.check)
        buttons.addStretch(1)
        buttons.addWidget(self.save)
        box.addLayout(buttons)
        outer.addWidget(self.box)
        self.notice = hint()
        self.notice.setVisible(False)
        outer.addWidget(self.notice)

        # -- jobs -----------------------------------------------------------
        self.table = QtWidgets.QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Model", "Task", "Where", "Status"])
        hh = self.table.horizontalHeader()
        hh.setStretchLastSection(True)
        for c in (0, 1, 2):
            hh.setSectionResizeMode(c, QtWidgets.QHeaderView.ResizeToContents)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.itemSelectionChanged.connect(self.select_job)
        self.table.cellClicked.connect(self.select_job)
        outer.addWidget(self.table, 1)
        jb = QtWidgets.QHBoxLayout()
        self.cancel = QtWidgets.QPushButton("Cancel job")
        self.cancel.setToolTip("Cancel the selected running job (scancel on the cluster). Finished results are kept.")
        self.cancel.clicked.connect(self.cancel_job)
        self.clear = QtWidgets.QPushButton("Clear selected")
        self.clear.setToolTip("Remove the selected finished jobs from the list (Ctrl- or Shift-click to pick "
                              "several). Records move to fem/out/jobs/cleared/; results are not touched.")
        self.clear.clicked.connect(self.clear_selected)
        jb.addWidget(self.cancel)
        jb.addWidget(self.clear)
        jb.addStretch(1)
        self.copy = QtWidgets.QToolButton()
        self.copy.setText("Copy record")
        self.copy.setToolTip("Copy the selected job's full record (JSON) to the clipboard")
        self.copy.clicked.connect(self.copy_record)
        jb.addWidget(self.copy)
        outer.addLayout(jb)

        # -- selected job ---------------------------------------------------
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        outer.addWidget(self.summary)
        self.details = QtWidgets.QPlainTextEdit()
        self.details.setReadOnly(True)
        self.details.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.details.setMaximumBlockCount(400)
        outer.addWidget(self.details, 1)

        self.load_profile()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(lambda: self.refresh() if self.isVisible() else None)
        self.timer.start(1000)
        self.refresh()

    def say(self, text):
        """The one notice line; hidden when there is nothing to say."""
        self.notice.setText(text)
        self.notice.setVisible(bool(text))

    # -- profiles -------------------------------------------------------------
    def current_key(self):
        return self.keys[self.run_on.index()]

    def use_profile(self, i):
        import FreeCAD as App
        self.load_profile()
        if App.ActiveDocument:
            try:
                self.service.select_profile(App.ActiveDocument, self.keys[i])
                self.say("")
            except ValueError as exc:
                self.say(str(exc))
        self.refresh()

    def load_profile(self):
        p = self.service.store.profiles[self.current_key()]
        remote = p["backend"] == "slurm"
        self.box.setTitle(p["name"])
        # show/hide, not a QStackedWidget: that sizes to its tallest page
        self.local_page.setVisible(not remote)
        self.remote_page.setVisible(remote)
        self.check.setVisible(remote)
        if remote:
            try:
                p = J.validate_profile(p)       # converts the old account/partition/CPU fields
            except ValueError as exc:
                self.say(str(exc))
            self.host.setText(p["host"])
            self.sbatch.setText(p.get("sbatch_options", ""))
        else:
            for w, v in ((self.mesh_threads, p["mesh_threads"]), (self.cpus, p["solve_cpus"]),
                         (self.memory, p["memory_gb"])):
                w.blockSignals(True)
                w.setValue(v)
                w.blockSignals(False)
        self.save.setEnabled(False)

    def mark_edited(self, *_):
        self.save.setEnabled(True)

    def edited_profile(self):
        p = dict(self.service.store.profiles[self.current_key()])
        if p["backend"] == "local":
            p.update(mesh_threads=self.mesh_threads.value(), solve_cpus=self.cpus.value(),
                     memory_gb=self.memory.value())
        else:
            p.update(host=self.host.text().strip(), sbatch_options=self.sbatch.text().strip())
        return J.validate_profile(p)

    def save_profile(self):
        try:
            self.service.store.save_profile(self.current_key(), self.edited_profile())
            self.save.setEnabled(False)
            self.say("Saved. Running jobs keep their own settings.")
        except (ValueError, OSError) as exc:
            self.say(str(exc))
        self.refresh()

    # -- jobs -----------------------------------------------------------------
    def selected_id(self):
        row = self.table.currentRow()
        item = self.table.item(row, 0) if row >= 0 else None
        return item.data(QtCore.Qt.UserRole) if item else None

    def selected_ids(self):
        return {self.table.item(i.row(), 0).data(QtCore.Qt.UserRole)
                for i in self.table.selectionModel().selectedRows()}

    def clearable(self):
        return {i for i in self.selected_ids()
                if self.service.records[i]["state"] in J.TERMINAL and i not in self.service.runners}

    def status(self, r, now):
        """(text, colour) for one job: one status, no separate percentage."""
        state, msg = r["state"], (r.get("message") or "")
        if state == "running":
            if msg.startswith("queued"):          # e.g. "queued on case-hpc (Priority)"
                reason = msg[msg.find("("):] if "(" in msg else ""
                return ("Queued " + reason).strip(), STATE_COLOUR["queued"]
            return (self.clock.text(r["id"], r["task"], int(r.get("progress") or 0), now,
                                    scs_stages.downloading(msg)), STATE_COLOUR["running"])
        if state == "queued":
            return "Queued", STATE_COLOUR["queued"]
        if state == "completed":
            return "Done · " + scs_stages.duration(r["updated"] - r["created"]), STATE_COLOUR["completed"]
        return state.capitalize(), STATE_COLOUR.get(state, "#666666")

    def refresh(self):
        import FreeCAD as App
        doc = App.ActiveDocument
        dock = self.parentWidget()
        if isinstance(dock, QtWidgets.QDockWidget):
            dock.setWindowTitle("SCS Job Manager — " + (doc.Label if doc else "no model open"))
        if doc is not None:
            key = self.service.profile_key(doc)
            if key in self.keys and self.keys[self.run_on.index()] != key and not self.save.isEnabled():
                self.run_on.set_index(self.keys.index(key))
                self.load_profile()
        self.run_on.setEnabled(doc is not None)
        now = time.time()
        selected = self.selected_id()
        chosen = self.selected_ids() & set(self.service.records) or ({selected} if selected else set())
        records = sorted(self.service.records.values(), key=lambda r: r["created"], reverse=True)
        self.table.blockSignals(True)
        self.table.setRowCount(len(records))
        rows = QtCore.QItemSelection()
        for row, r in enumerate(records):
            text, colour = self.status(r, now)
            job = (r.get("remote") or {}).get("job_id")
            where = r["profile"]["name"].replace(" workstation", "") + (" · job %s" % job if job else "")
            for col, value in enumerate((r["document"]["label"], r["task"], where, text)):
                item = QtWidgets.QTableWidgetItem(value)
                item.setData(QtCore.Qt.UserRole, r["id"])
                if col == 3:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(colour)))
                    item.setToolTip(r.get("message") or r["state"])
                self.table.setItem(row, col, item)
            # Rows move when a job is added, so the selection follows job ids.
            if r["id"] in chosen:
                rows.select(self.table.model().index(row, 0), self.table.model().index(row, 3))
            if selected == r["id"]:
                self.table.selectionModel().setCurrentIndex(self.table.model().index(row, 0),
                                                            QtCore.QItemSelectionModel.NoUpdate)
        self.table.selectionModel().select(rows, QtCore.QItemSelectionModel.ClearAndSelect)
        if selected is None and records:
            self.table.selectRow(0)
        self.table.blockSignals(False)
        self.show_details()

    def select_job(self, *_):
        self.probe_text = None
        self.show_details()

    def set_details(self, text):
        # The 1 s refresh must not reset the scroll position or cursor.
        if text == self.details.toPlainText():
            return
        bar = self.details.verticalScrollBar()
        at_end = bar.value() == bar.maximum()
        pos = bar.value()
        self.details.setPlainText(text)
        bar.setValue(bar.maximum() if at_end and pos else pos)

    def show_details(self):
        ident = self.selected_id()
        runner = self.service.runners.get(ident)
        self.cancel.setEnabled(bool(runner and runner.is_running()))
        self.clear.setEnabled(bool(self.clearable()))
        self.copy.setEnabled(ident is not None)
        if self.probe_text is not None:
            self.summary.setText("")
            return self.set_details(self.probe_text)
        if not ident:
            self.summary.setText("")
            return self.set_details("")
        r = self.service.records[ident]
        p = r["profile"]
        where = ("%s, sbatch %s" % (p.get("host"), p.get("sbatch_options"))) if p.get("backend") == "slurm" \
            else p.get("name", "local")
        rows = [("Model", r["document"]["label"]),
                ("Run", (r.get("args") or ["-"])[0]),
                ("On", where),
                ("Submitted", time.strftime("%a %d %b %H:%M", time.localtime(r["created"]))
                 + " · " + scs_stages.duration((time.time() if r["state"] not in J.TERMINAL else r["updated"])
                                                    - r["created"])),
                ("State", r.get("message") or r["state"])]
        if r["state"] in ("running", "cancelling"):
            rows.append(("Resources", self.monitor.line(ident) or "sampling\u2026"))
        self.summary.setText("<table>" + "".join(
            "<tr><td style='color:gray; padding-right:8px'>%s</td><td>%s</td></tr>" % (k, v) for k, v in rows)
            + "</table>")
        path = r.get("log_path")
        tail = ""
        if path and Path(path).is_file():
            with open(path, "rb") as f:
                f.seek(max(0, Path(path).stat().st_size - 12000))
                tail = f.read().decode("utf8", "replace")
                lines = tail.splitlines()
                if f.tell() > 12000:
                    lines = lines[1:]              # the first line is cut by the seek
                # the size field prints one line per chunk; keep the log readable
                lines = [l for l in lines if not l.startswith(("PROGRESS ", "  field "))]
                tail = "\n".join(lines[-150:])
        self.set_details(tail)

    def copy_record(self):
        ident = self.selected_id()
        if ident:
            QtWidgets.QApplication.clipboard().setText(json.dumps(self.service.records[ident], indent=2))
            self.say("Record copied.")

    def cancel_job(self):
        ident = self.selected_id()
        if ident:
            self.service.cancel(ident)
        self.refresh()

    def clear_selected(self):
        idents = self.clearable()
        self.table.clearSelection()
        self.details.clear()
        count = self.service.clear_finished(idents)
        self.say("Cleared %d job(s); records are in %s. Results were not touched."
                 % (count, self.service.store.root/"cleared"))
        self.refresh()

    # -- cluster check (read only) -------------------------------------------
    def check_connection(self):
        try:
            self.probe_profile = self.edited_profile()
            args = J.ssh_probe_args(self.probe_profile)
        except ValueError as exc:
            self.say(str(exc))
            return
        self.probe = QtCore.QProcess(self)
        self.probe.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.probe.finished.connect(self.probe_finished)
        self.probe.errorOccurred.connect(self.probe_error)
        self.check.setEnabled(False)
        self.say("Checking the cluster…")
        self.probe.start("ssh", args)
        # The cluster Python's first import of numpy/scipy/gmsh/vtk from the
        # NFS home takes ~16 s when cold (5 s warm), measured 2026-10-15; a 20 s
        # limit killed healthy checks and reported "exit code 9".
        self.probe_timed_out = False
        self.probe_timeout = QtCore.QTimer(self)
        self.probe_timeout.setSingleShot(True)
        self.probe_timeout.timeout.connect(self.probe_expired)
        self.probe_timeout.start(PROBE_TIMEOUT_MS)
        self.probe_text = "Checking the cluster…"
        self.show_details()

    def probe_expired(self):
        self.probe_timed_out = True
        self.probe.kill()

    def probe_error(self, error):
        self.check.setEnabled(True)
        self.say("SSH check: " + self.probe.errorString())
        self.probe_text = "SSH check failed: " + self.probe.errorString()
        self.show_details()

    def probe_finished(self, code, status):
        self.probe_timeout.stop()
        self.check.setEnabled(True)
        output = bytes(self.probe.readAllStandardOutput()).decode("utf8", "replace")
        verdict, details = J.describe_probe(code, output, self.probe_profile)
        if self.probe_timed_out:
            verdict = ("Timed out after %d s, with the steps above done. The usual cause is the cluster Python "
                       "importing numpy/scipy/gmsh/vtk from a cold NFS cache; check again." % (PROBE_TIMEOUT_MS//1000))
        self.say(verdict)
        self.probe_text = verdict + "\n\n" + details + "\n\nNothing was submitted. Click a job to return to it."
        self.show_details()

    def shutdown(self):
        self.timer.stop()


def show():
    """Open the Job Manager. Re-running rebuilds the window from the current
    code (the jobs live in scs_job_service, which is never rebuilt)."""
    import FreeCADGui as Gui
    mw = Gui.getMainWindow()
    service().adopt_open_panels()
    dock = mw.findChild(QtWidgets.QDockWidget, "SCS_JobManagerDock")
    if dock is None:
        dock = QtWidgets.QDockWidget("SCS Job Manager", mw)
        dock.setObjectName("SCS_JobManagerDock")
        mw.addDockWidget(QtCore.Qt.RightDockWidgetArea, dock)
    old = dock.widget()
    if old is not None and hasattr(old, "timer"):
        old.timer.stop()
    dock.setWidget(JobManager())
    if old is not None:
        old.deleteLater()
    dock.show()
    dock.raise_()
    return dock
