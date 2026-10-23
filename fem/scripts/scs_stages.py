"""Job progress as named stages, for the Mesh Generator and the Job Manager.

The mesh worker reports PROGRESS 0-100, but the number is not time: 0-60 is
the size field (real progress, chunk by chunk), 60->75 one jump when gmsh's
tetrahedralization ends, then tissues 75-85, the dura check 85-92 and the
preview 92-100. A 415 M-tet build sat at 60 % for gmsh and at 75 % for the
tissues, and its 4.9 h time estimate (a workstation calibration) was wrong by
20x for gmsh on 24 smp cores (188 s), so neither the percentage nor a clock
against that estimate says where a job is (Mohamed, 2026-10-23: "I'm not sure
why it's reading 28 %, I only started it about 5 minutes ago"). This shows
the stage, how long it has been running, a bar only where the worker reports
real progress (busy otherwise), and an ETA only where it can be measured.

UI only: nothing here is part of the mesh or solve signature.
"""
import json
import os
import time
from pathlib import Path

try:
    from PySide import QtCore, QtWidgets
except ImportError:                                     # pragma: no cover
    from PySide6 import QtCore, QtWidgets

# (start %, name, reports progress inside the stage) -- mesh_preview.py's scale
MESH_STAGES = ((0, "Size field", True), (60, "Tetrahedra", False), (75, "Tissues", False),
               (85, "Dura check", False), (92, "Preview", False))
# solve_lead.py's scale: tissues and islands to 5, the Elmer mesh to 12, then
# one step per contact to 95; saving (and on a cluster the surface basis) after.
SOLVE_STAGES = ((0, "Tissues", False), (5, "Elmer mesh", False), (12, "Contacts", True),
                (95, "Save", False))
STAGES = {"Mesh": MESH_STAGES, "Solve": SOLVE_STAGES}
DOWNLOAD = "Download"      # a cluster job's results coming back, after its last stage


def downloading(message):
    """True while the job service is copying a cluster job's results here."""
    return "fetching results" in (message or "")

DONE_FG, DONE_BG = "#1b5e20", "#c8e6c9"
RUN_FG, RUN_BG = "#7a4f01", "#ffe0b2"
FAIL_FG, FAIL_BG = "#8e0000", "#ffcdd2"
IDLE_FG = "#8a8a8a"        # stages not reached yet: readable on light and dark themes


def duration(s):
    s = max(0, int(s))
    if s < 90:
        return "%d s" % s
    if s < 3600:
        return "%d min" % round(s/60)
    return "%d h %02d min" % (s//3600, (s % 3600)//60)


def stage_of(task, pct, download=False):
    """(index, name, fraction within the stage or None, stages) for a task."""
    stages = STAGES.get(task)
    if not stages:
        return None, None, max(0, min(1, pct/100.0)), None
    if download:
        return len(stages), DOWNLOAD, None, stages
    k = max(i for i, (start, _n, _r) in enumerate(stages) if pct >= start)
    start, name, reports = stages[k]
    end = stages[k + 1][0] if k + 1 < len(stages) else 100
    return k, name, ((pct - start)/(end - start) if reports and end > start else None), stages


STAGE_FILE = Path(__file__).resolve().parents[1] / "out" / "jobs" / "stages.json"


class StageClock:
    """When each job entered its current stage, and its rate there (for an
    ETA where the stage reports progress). Saved to STAGE_FILE whenever a job
    changes stage, so closing and reopening a panel, or FreeCAD, keeps the
    times (Mohamed: the status time "resets to zero when I close and open the
    macro window"). The Job Manager and the Mesh Generator share the file. A
    stage that began while no panel was watching counts from when one first
    saw it."""

    def __init__(self, path=STAGE_FILE):
        self.path = Path(path) if path else None
        self.jobs = {}
        if self.path and self.path.is_file():
            try:
                self.jobs = json.loads(self.path.read_text())
            except (OSError, ValueError):
                self.jobs = {}

    def _save(self):
        if not self.path:
            return
        try:
            # merge with what another panel may have written since
            disk = json.loads(self.path.read_text()) if self.path.is_file() else {}
            disk.update(self.jobs)
            old = time.time() - 30*86400
            self.jobs = {k: j for k, j in disk.items() if j.get("since", 0) > old}
            tmp = self.path.with_suffix(".tmp")
            tmp.write_text(json.dumps(self.jobs))
            os.replace(tmp, self.path)
        except (OSError, ValueError):
            pass

    def update(self, ident, task, pct, now=None, download=False):
        now = now or time.time()
        k, name, frac, _ = stage_of(task, pct, download)
        j = self.jobs.get(ident)
        if j is None or j["stage"] != k:
            if self.path and self.path.is_file():       # another panel may know better
                try:
                    j = json.loads(self.path.read_text()).get(ident) or j
                except (OSError, ValueError):
                    pass
            if j is None or j["stage"] != k:
                j = dict(stage=k, since=now, f0=frac, t0=now)
                self.jobs[ident] = j
                self._save()
            else:
                self.jobs[ident] = j
        j["frac"], j["now"] = frac, now
        return j

    def eta(self, ident):
        j = self.jobs.get(ident)
        if not j or j.get("frac") is None or j.get("f0") is None or j["now"] - j["t0"] < 20:
            return None
        df = j["frac"] - j["f0"]
        if df <= 0.02:
            return None
        return (1 - j["frac"])*(j["now"] - j["t0"])/df

    def text(self, ident, task, pct, now=None, download=False):
        """'Tissues · 4 min' or 'Size field · 40 % · ~3 min left'."""
        now = now or time.time()
        j = self.update(ident, task, pct, now, download)
        k, name, frac, _ = stage_of(task, pct, download)
        if name is None:
            return "%d %%" % pct
        parts = [name, duration(now - j["since"])]
        if frac is not None:
            parts.insert(1, "%d %%" % (100*frac))
            eta = self.eta(ident)
            if eta is not None:
                parts.append("~%s left" % duration(eta))
        return " \u00b7 ".join(parts)


class StageView(QtWidgets.QWidget):
    """Stage markers, one thin bar, one caption line."""

    def __init__(self, task="Mesh", parent=None):
        super().__init__(parent)
        self.task, self.clock, self.ident = task, StageClock(None), "panel"
        self._pct, self._started, self._dl = 0, None, False
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(3)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(2)
        self.marks = []
        self.names = [name for _s, name, _r in STAGES.get(task, ())]
        if self.names:
            self.names.append(DOWNLOAD)                 # shown only once a cluster job's results come back
        for name in self.names:
            m = QtWidgets.QLabel(name)
            m.setAlignment(QtCore.Qt.AlignCenter)
            m.setMinimumHeight(22)
            row.addWidget(m, 1)
            self.marks.append(m)
        lay.addLayout(row)
        self.bar = QtWidgets.QProgressBar()
        self.bar.setTextVisible(False)
        self.bar.setMaximumHeight(6)
        lay.addWidget(self.bar)
        self.caption = QtWidgets.QLabel("")
        self.caption.setWordWrap(True)
        lay.addWidget(self.caption)
        self.resources = QtWidgets.QLabel("")          # scs_monitor's line, while running
        self.resources.setStyleSheet("color: gray;")
        self.resources.setToolTip("Live resource use of this job (scs_monitor): resident memory now and at "
                                  "peak against the allocation, and how many of its cores are busy, from "
                                  "the job's cgroup on its cluster node or /proc here; refreshed every second.")
        lay.addWidget(self.resources)
        self.idle()

    def _paint(self, current, state):
        idle_fg = IDLE_FG
        for i, m in enumerate(self.marks):
            if state == "done" or (current is not None and i < current):
                fg, bg, bold = DONE_FG, DONE_BG, False
            elif current is not None and i == current:
                fg, bg, bold = (FAIL_FG, FAIL_BG, True) if state == "failed" else (RUN_FG, RUN_BG, True)
            else:
                fg, bg, bold = idle_fg, "transparent", False
            m.setStyleSheet("QLabel { color: %s; background: %s; border-radius: 3px; padding: 1px 4px; %s }"
                            % (fg, bg, "font-weight: bold;" if bold else ""))

    def _download_mark(self, visible):
        self._dl = visible
        if self.marks:
            self.marks[-1].setVisible(visible)

    def idle(self, text=""):
        self._started = None
        self._download_mark(False)
        self.resources.setVisible(False)
        self._paint(None, "idle")
        self.bar.setRange(0, 1)
        self.bar.setValue(0)
        self.bar.setVisible(False)
        self.set_caption(text)

    def attach(self, ident, submitted):
        """Follow a job: its stage times come from the shared, saved clock and
        the total from when it was submitted, so reopening the panel keeps both."""
        if ident and ident != self.ident:
            self.ident, self.clock = ident, StageClock()
        if submitted:
            self._started = submitted

    def start(self, text=""):
        self._started = time.time()
        self.clock, self.ident = StageClock(None), "panel"
        self._download_mark(False)
        self.bar.setVisible(True)
        self.progress(0)
        self.set_caption(text)

    def progress(self, pct):
        self._pct = pct
        k, _name, frac, _ = stage_of(self.task, pct)
        self.clock.update(self.ident, self.task, pct)
        self._paint(k, "running")
        if frac is None:
            self.bar.setRange(0, 0)               # busy: the worker reports nothing here
        else:
            self.bar.setRange(0, 1000)
            self.bar.setValue(int(1000*frac))

    def tick(self, note=""):
        """Call every second or so while running: refreshes stage time / ETA."""
        if self._started is None:
            return
        dl = downloading(note) or self._dl             # once downloading, until done/failed
        text = self.clock.text(self.ident, self.task, self._pct, download=dl)    # "Tissues · 4 min"
        if self.marks:
            k = stage_of(self.task, self._pct, dl)[0]
            if dl:
                self._download_mark(True)
                self._paint(k, "running")
                self.bar.setRange(0, 0)
            for i, n in enumerate(self.names):
                self.marks[i].setText(text if i == k else n)
        total = duration(time.time() - self._started)
        note = "" if dl else note                       # the Download mark already says it
        self.set_caption(("%s · %s total" % (note, total)) if note else "%s total" % total)

    def done(self, text):
        self._started = None
        self.resources.setVisible(False)
        for i, n in enumerate(self.names):
            self.marks[i].setText(n)
        self._paint(None, "done")
        self.bar.setVisible(False)
        self.set_caption(text)

    def failed(self, text):
        self._started = None
        self.resources.setVisible(False)
        k, _n, _f, _ = stage_of(self.task, self._pct, self._dl)
        self._paint(k, "failed")
        self.bar.setVisible(False)
        self.set_caption(text, FAIL_FG)

    def set_resources(self, text):
        self.resources.setText(text)
        self.resources.setVisible(bool(text) and self._started is not None)

    def set_caption(self, text, colour=None):
        self.caption.setText(text)
        self.caption.setVisible(bool(text))
        self.caption.setStyleSheet("color: %s;" % (colour or "gray"))
