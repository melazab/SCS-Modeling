"""Shared Qt subprocess runner with progress, persistent logs and cancellation.

Scripts emit PROGRESS <0-100>, and one final JSON report with the specified
prefix. Full stdout/stderr and completion status are kept alongside artifacts.
Silence produces a warning, not an automatic kill: sorting/AMG may be slow.
"""
import json
import os
import time

try:
    from PySide import QtCore
except ImportError:                                     # pragma: no cover
    try:
        from PySide6 import QtCore
    except ImportError:
        from PySide2 import QtCore

# Long native phases can legitimately be silent; never infer a deadlock from
# stdout timing alone. Record a diagnostic after 15 minutes, leaving Abort
# available. Actual crashes/nonzero exits still fail immediately.
STALL_TIMEOUT_MS = 900_000

STATUS_COLORS = {
    "idle":    ("#333333", "#e0e0e0"),
    "running": ("#7a4f01", "#ffe0b2"),
    "done":    ("#1b5e20", "#c8e6c9"),
    "failed":  ("#8e0000", "#ffcdd2"),
    "aborted": ("#7a4f01", "#ffe0b2"),
    "gated":   ("#7a4f01", "#ffe0b2"),
    "active":  ("#1b5e20", "#c8e6c9"),
}


def set_status(label, level, text):
    fg, bg = STATUS_COLORS.get(level, ("#333333", "#e0e0e0"))
    label.setText(text)
    label.setStyleSheet(
        "background: %s; color: %s; border-radius: 3px; font-weight: bold;" % (bg, fg))


class AsyncProcessRunner(QtCore.QObject):
    """Runs one `python_exe script <args>` subprocess at a time, with a
    watchdog, a verified-dead Abort, an elapsed-time tick, and LIVE
    progress parsed from the subprocess's own stdout.

    Signals:
        tick(float)       -- elapsed seconds, once a second while running
        progress(int)     -- 0-100, parsed from "PROGRESS <int>" lines,
                              emitted incrementally as they arrive
        finished_ok(dict) -- the parsed JSON report line, plus
                              "_stdout_tail" (last 4000 chars of stdout)
        failed(str)       -- message
        aborted(str)      -- message
    """
    tick = QtCore.Signal(float)
    progress = QtCore.Signal(int)
    stage_message = QtCore.Signal(str)
    finished_ok = QtCore.Signal(dict)
    failed = QtCore.Signal(str)
    aborted = QtCore.Signal(str)

    def __init__(self, parent, python_exe, script, json_prefix,
                 stall_timeout_ms=STALL_TIMEOUT_MS):
        super(AsyncProcessRunner, self).__init__(parent)
        self.python_exe = python_exe
        self.script = script
        self.json_prefix = json_prefix
        self.stall_timeout_ms = stall_timeout_ms
        self._last_output = None
        self._stalled = False
        self.process = None
        self.start_time = None
        self._run_token = 0
        self._aborted = False
        self._stdout_buf = ""
        self._line_buf = ""
        self._last_progress = -1
        self._stderr_buf = ""
        self.log_path = None
        self.delegate = None   # a RemoteJob driving this runner's signals, see attach()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._on_tick)

    def is_running(self):
        if self.delegate is not None:
            return self.delegate.is_running()
        return self.process is not None and self.process.state() != QtCore.QProcess.NotRunning

    def attach(self, job):
        """Let `job` (scs_remote.RemoteJob: same signals, runs elsewhere) drive
        this runner, so a panel's own handlers, busy state and Abort button
        work unchanged for a job that runs on a cluster."""
        self.delegate = job
        self._last_progress = -1
        self.log_path = job.log_path
        for name in ('tick', 'progress', 'stage_message'):
            getattr(job, name).connect(getattr(self, name).emit)
        job.progress.connect(lambda pct: setattr(self, '_last_progress', pct))

        def relay(signal):
            def emit(value):
                if self.delegate is job:
                    self.delegate = None
                signal.emit(value)
            return emit
        job.finished_ok.connect(relay(self.finished_ok))
        job.failed.connect(relay(self.failed))
        job.aborted.connect(relay(self.aborted))

    def start(self, args, env_overrides=None):
        if self.is_running():
            return
        self.delegate = None
        proc = QtCore.QProcess(self)
        env = QtCore.QProcessEnvironment.systemEnvironment()
        for k, v in (env_overrides or {}).items():
            env.insert(k, str(v))
        env.insert("PYTHONUNBUFFERED", "1")
        proc.setProcessEnvironment(env)
        proc.setProgram(self.python_exe)
        proc.setArguments([self.script] + list(args))
        self._aborted = False
        self._stdout_buf = ""
        self._stderr_buf = ""
        self._line_buf = ""
        self._last_progress = -1
        log_dir = str(args[0]) if args and os.path.isdir(str(args[0])) else os.path.dirname(self.script)
        self.log_path = os.path.join(log_dir, os.path.basename(self.script) + '.log')
        with open(self.log_path, 'w') as f:
            f.write('Started ' + time.strftime('%Y-%m-%d %H:%M:%S') + '\n')
        self._run_token += 1
        token = self._run_token
        proc.readyReadStandardOutput.connect(lambda t=token: self._on_stdout(t))
        proc.readyReadStandardError.connect(lambda t=token: self._on_stderr(t))
        proc.finished.connect(lambda code, status, t=token: self._on_finished(code, status, t))
        proc.errorOccurred.connect(lambda err, t=token: self._on_error(err, t))
        self.process = proc
        self.start_time = time.time()
        self._last_output = time.time()
        self._stalled = False
        proc.start()
        self.timer.start(1000)

    def abort(self):
        if not self.is_running():
            return
        if self.delegate is not None:
            return self.delegate.abort()
        self._aborted = True
        pid = self.process.processId()
        if pid:
            try:
                import psutil
                for child in reversed(psutil.Process(pid).children(recursive=True)):
                    try:
                        child.kill()
                    except psutil.NoSuchProcess:
                        pass
            except ImportError:
                pass
            except psutil.NoSuchProcess:
                pass
        # Local Elmer workers create a process group containing their MPI ranks.
        if pid and hasattr(os, 'killpg'):
            try:
                if os.getpgid(pid) == pid:
                    import signal
                    os.killpg(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        self.process.kill()   # SIGKILL on Unix -- immediate, uncatchable
        self.process.waitForFinished(5000)
        alive = self._pid_alive(pid)
        self.timer.stop()
        self.start_time = None
        why = "aborted by user"
        self._log("\n" + why + "\n")
        if alive:
            self.aborted.emit("%s, but pid %d is still alive 5s after SIGKILL -- "
                              "investigate by hand (ps -p %d)" % (why, pid, pid))
        else:
            self.aborted.emit("%s -- pid %d confirmed dead" % (why, pid))

    @staticmethod
    def _pid_alive(pid):
        if not pid:
            return False
        try:
            import psutil
            return psutil.pid_exists(pid)
        except ImportError:
            try:
                os.kill(pid, 0)
                return True
            except ProcessLookupError:
                return False
            except PermissionError:
                return True

    def _on_tick(self):
        if self.start_time is None:
            return
        now = time.time()
        self.tick.emit(now - self.start_time)
        if (self._last_output is not None
                and (now - self._last_output) * 1000.0 > self.stall_timeout_ms
                and self.is_running()):
            if not self._stalled:
                self._stalled = True
                self._log('\nNo recent output; worker still running. Use Abort to cancel.\n')

    def _log(self, chunk):
        if self.log_path:
            with open(self.log_path, 'a') as f:
                f.write(chunk)

    def _on_stderr(self, token):
        if token != self._run_token or self.process is None:
            return
        chunk = bytes(self.process.readAllStandardError()).decode('utf-8', 'replace')
        if chunk:
            self._last_output = time.time()
            self._stderr_buf = (self._stderr_buf + chunk)[-16000:]
            self._log(chunk)

    def _on_stdout(self, token):
        if token != self._run_token or self.process is None:
            return
        chunk = bytes(self.process.readAllStandardOutput()).decode("utf-8", "replace")
        # Any output at all resets the stall clock -- this is what distinguishes
        # "slow but working" from "hung", see STALL_TIMEOUT_MS above.
        self._last_output = time.time()
        self._log(chunk)
        self._stdout_buf = (self._stdout_buf + chunk)[-64000:]
        self._line_buf += chunk
        while "\n" in self._line_buf:
            line, self._line_buf = self._line_buf.split("\n", 1)
            if line.startswith(('RESOURCE ', 'Contact ', 'Active:', '  classified ')):
                self.stage_message.emit(line.strip())
            if line.startswith("PROGRESS "):
                try:
                    pct = int(line.split(None, 1)[1].strip())
                except (IndexError, ValueError):
                    continue
                pct = max(0, min(100, pct))
                if pct != self._last_progress:
                    self._last_progress = pct
                    self.progress.emit(pct)


    def _on_error(self, error, token):
        if token != self._run_token:
            return
        if not self.is_running():
            self.timer.stop()
            self.start_time = None
            if not self._aborted and error == QtCore.QProcess.FailedToStart:
                self.failed.emit("Could not start worker: %s. Log: %s" % (self.process.errorString(), self.log_path))

    def _on_finished(self, exit_code, exit_status, token):
        if token != self._run_token:
            return
        self.timer.stop()
        self.start_time = None
        self._on_stdout(token)
        self._on_stderr(token)
        stderr = self._stderr_buf
        self._log('\nExit %s, status %s\n' % (exit_code, exit_status))
        stdout = self._stdout_buf
        if self._aborted:
            return  # abort() already emitted `aborted`
        if exit_status == QtCore.QProcess.CrashExit or exit_code != 0:
            reason = ('worker killed (SIGKILL); check the kernel OOM log' if exit_code in (9, -9, 137)
                      else 'worker failed (exit %s)' % exit_code)
            self.failed.emit("%s\nLog: %s\n%s" % (reason, self.log_path, stderr[-4000:] or stdout[-4000:]))
            return
        line = next((ln for ln in stdout.splitlines() if ln.startswith(self.json_prefix)), None)
        if line is None:
            self.failed.emit("no %s in output\n%s" % (self.json_prefix.strip(), stdout[-4000:]))
            return
        try:
            report = json.loads(line[len(self.json_prefix):])
            if not isinstance(report, dict):
                raise ValueError('report must be an object')
        except ValueError as exc:
            self.failed.emit('Invalid worker report: %s. Log: %s' % (exc, self.log_path))
            return
        report["_stdout_tail"] = stdout[-4000:]
        report["log_path"] = self.log_path
        self.progress.emit(100)
        self.finished_ok.emit(report)
