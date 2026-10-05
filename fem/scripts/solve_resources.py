"""Local worker budget, process-tree monitoring and stage memory records."""
import os
import math
import threading
import psutil


def available_gb():
    return psutil.virtual_memory().available / 1e9


class Resources:
    def __init__(self, cpus=None, memory_gb=None):
        available_cpus = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
        self.cpus = int(cpus if cpus is not None else os.environ.get('OMP_NUM_THREADS', '1'))
        if not 1 <= self.cpus <= available_cpus:
            raise ValueError('CPUs must be between 1 and %d' % available_cpus)
        available = available_gb()
        requested = float(memory_gb if memory_gb is not None else os.environ.get('SCS_SOLVE_MEMORY_GB', '0'))
        if not math.isfinite(requested) or requested < 0:
            raise ValueError('Memory budget must be nonnegative (0 = automatic)')
        self.memory_gb = min(requested, .8 * available) if requested else .8 * available
        self.peak_gb = 0.
        self.stage = 'loading'
        self.records = []
        self.stop = threading.Event()
        self.process = psutil.Process()

    def usage(self):
        # Sum RSS conservatively, including MPI launchers and their workers.
        processes = [self.process] + self.process.children(recursive=True)
        total = 0
        for process in processes:
            try:
                total += process.memory_info().rss
            except psutil.NoSuchProcess:
                pass
        return total / 1e9

    def mark(self, stage):
        usage = self.usage()
        self.peak_gb = max(self.peak_gb, usage)
        self.records.append(dict(stage=stage, rss_gb=round(usage, 3), peak_gb=round(self.peak_gb, 3)))
        self.stage = stage
        print('RESOURCE %s: %.3f GB RSS, %.3f GB peak, %.3f GB budget' %
              (stage, usage, self.peak_gb, self.memory_gb), flush=True)

    def _watch(self):
        while not self.stop.wait(.25):
            usage = self.usage()
            self.peak_gb = max(self.peak_gb, usage)
            if usage > self.memory_gb:
                print('Memory budget exceeded during %s: %.2f GB > %.2f GB. '
                      'Stopping worker and child solvers; no solution published.' %
                      (self.stage, usage, self.memory_gb), flush=True)
                for child in reversed(self.process.children(recursive=True)):
                    try:
                        child.kill()
                    except psutil.NoSuchProcess:
                        pass
                os._exit(2)

    def __enter__(self):
        os.environ['SCS_SOLVE_MEMORY_GB'] = str(self.memory_gb)
        self.thread = threading.Thread(target=self._watch, daemon=True)
        self.thread.start()
        self.mark('worker start')
        return self

    def __exit__(self, *exc):
        self.stop.set()
        self.thread.join(timeout=1)
