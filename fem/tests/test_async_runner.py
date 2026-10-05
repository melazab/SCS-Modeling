"""Real subprocess lifecycle tests; no FreeCAD document is needed."""
from pathlib import Path
import sys
import tempfile
import time
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    from async_runner import AsyncProcessRunner, QtCore
except ImportError:
    AsyncProcessRunner = QtCore = None


@unittest.skipIf(QtCore is None, 'Qt subprocess tests require PySide from FreeCAD or PySide6')
class RunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QtCore.QCoreApplication.instance() or QtCore.QCoreApplication([])

    def run_worker(self, code, silence=False, abort=False):
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / 'worker.py'
            script.write_text(code)
            runner = AsyncProcessRunner(None, sys.executable, str(script), 'RESULT ', stall_timeout_ms=10)
            loop = QtCore.QEventLoop()
            outcomes, progress = [], []
            runner.progress.connect(progress.append)
            def finish(kind, data):
                outcomes.append((kind, data))
                loop.quit()
            runner.finished_ok.connect(lambda data: finish('ok', data))
            runner.failed.connect(lambda data: finish('failed', data))
            runner.aborted.connect(lambda data: finish('aborted', data))
            runner.start([tmp])
            if abort:
                QtCore.QTimer.singleShot(300, runner.abort)
            if silence:
                def tick():
                    runner._last_output = time.time() - 1
                    runner._on_tick()
                QtCore.QTimer.singleShot(25, tick)
            timeout = QtCore.QTimer()
            timeout.setSingleShot(True)
            timeout.timeout.connect(loop.quit)
            timeout.start(5000)
            loop.exec()
            timeout.stop()
            if runner.is_running():
                runner.abort()
                self.fail('Worker did not finish in five seconds')
            self.assertEqual(len(outcomes), 1)
            log = Path(runner.log_path).read_text()
            return outcomes[0], progress, log

    def test_drains_final_output_and_keeps_stderr(self):
        outcome, progress, log = self.run_worker(
            "import sys\nprint('PROGRESS 85')\nprint('diagnostic',file=sys.stderr)\n"
            "sys.stdout.write('RESULT {\"answer\": 42}')\n")
        self.assertEqual(outcome[0], 'ok')
        self.assertEqual(outcome[1]['answer'], 42)
        self.assertIn(85, progress)
        self.assertIn('diagnostic', log)
        self.assertIn('Exit 0', log)

    def test_nonzero_exit_is_not_success(self):
        outcome, _, log = self.run_worker("import sys\nprint('explanation',file=sys.stderr)\nsys.exit(7)\n")
        self.assertEqual(outcome[0], 'failed')
        self.assertIn('explanation', outcome[1])
        self.assertIn('Exit 7', log)

    def test_silence_does_not_kill_a_working_process(self):
        outcome, _, log = self.run_worker("import time\ntime.sleep(.15)\nprint('RESULT {}')", silence=True)
        self.assertEqual(outcome[0], 'ok')
        self.assertIn('worker still running', log)

    def test_bad_report_is_a_reported_failure(self):
        outcome, _, _ = self.run_worker("print('RESULT {invalid json}')")
        self.assertEqual(outcome[0], 'failed')
        self.assertIn('Invalid worker report', outcome[1])

    def test_abort_stops_solver_child(self):
        import psutil
        outcome, _, log = self.run_worker(
            "import os,subprocess,sys,time\nos.setsid()\n"
            "child=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'])\n"
            "print('CHILD_PID',child.pid,flush=True)\ntime.sleep(30)\n", abort=True)
        self.assertEqual(outcome[0], 'aborted')
        pid = int(next(line for line in log.splitlines() if line.startswith('CHILD_PID')).split()[1])
        self.assertTrue(not psutil.pid_exists(pid) or psutil.Process(pid).status() == psutil.STATUS_ZOMBIE)
