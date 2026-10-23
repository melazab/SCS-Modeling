"""Live resource use of running SCS jobs, every second: a lightweight htop for
the panels.

Mohamed (2026-10-23) asked for live memory and CPU in the Mesh Generator
(and later the Potential Visualizer), then for 1 s updates. The first reading
made the point: the 415 M-tet mesh spent ~20 min on one core (writing a 26 GB
mesh.msh nothing reads, then compressing mesh.npz) while 23 of its 24 cores
sat idle.

Cluster jobs: SLURM's own accounting (sstat) refreshes about every 30 s, and
an ssh plus a SLURM query every second would load the controller. Instead,
one ssh stays open to the job's node for as long as the job runs, and a small
shell loop there reads the job's cgroup (v1 on Pioneer) once a second:
memory.usage_in_bytes minus memory.stat's total_inactive_file (the working
set, the usual "memory in use"; forked workers' shared pages counted once),
memory.max_usage_in_bytes (peak), memory.limit_in_bytes (the allocation),
cpuacct.usage (CPU time, so busy cores = its change over wall time) and
cpuset.cpus (cores allocated). The node comes from sacct, asked once.

Local jobs: /proc for the worker and its descendants, every second (PSS from
smaps_rollup where the kernel has it, so forked workers are not counted
several times; VmRSS otherwise).

sample(ident) -> dict(rss_gb, peak_gb, alloc_gb, cores_busy, cores, node, t, age_s)
line(ident)   -> "RAM 85 / 600 GB (peak 130) · 1.0 of 24 cores busy · smpt12"

UI only: nothing here touches a job, its inputs or its signature.
"""
import os
import re
import time

try:
    from PySide import QtCore
except ImportError:                                     # pragma: no cover
    from PySide6 import QtCore

TICK_MS = 1000
NODE_RETRY_S = 15
CLK = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100

# Runs on the compute node: one line a second while the job's cgroup exists.
NODE_LOOP = r'''
u=$(id -u); j=%(job)s; g=/sys/fs/cgroup
m=$g/memory/slurm/uid_$u/job_$j; c=$g/cpuacct/slurm/uid_$u/job_$j
n=$(tr ',' '\n' < $g/cpuset/slurm/uid_$u/job_$j/cpuset.cpus 2>/dev/null | awk -F- '{n+=($2==""?1:$2-$1+1)} END{print n+0}')
while [ -d "$m" ]; do
  echo "S $(date +%%s.%%N) $(cat $m/memory.usage_in_bytes) $(awk '/^total_inactive_file /{print $2}' $m/memory.stat) $(cat $m/memory.max_usage_in_bytes) $(cat $m/memory.limit_in_bytes) $(cat $c/cpuacct.usage) $n"
  sleep 1
done
echo END
'''


def _gb(text):
    """'85193668K' / '600G' / '2048M' -> GB."""
    m = re.match(r"([\d.]+)([KMGT]?)", text or "")
    if not m:
        return None
    return float(m.group(1))*{"": 1e-9, "K": 1e-6, "M": 1e-3, "G": 1.0, "T": 1e3}[m.group(2)]


class Monitor(QtCore.QObject):
    updated = QtCore.Signal()

    def __init__(self):
        super().__init__()
        self.latest = {}          # ident -> sample dict
        self.nodes = {}           # ident -> (node, cores, alloc_gb)
        self.streams = {}         # ident -> QProcess streaming the node loop
        self._opened = {}         # ident -> when its stream was last opened
        self.finished = set()     # idents whose job has ended
        self._asked = {}          # ident -> last time sacct was asked
        self._cpu = {}            # local: ident -> (cpu_s, t)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.poll)
        self.timer.start(TICK_MS)
        QtCore.QTimer.singleShot(0, self.poll)

    def jobs(self):
        import scs_job_service
        return [r for r in scs_job_service.service().records.values() if r["state"] in ("running", "cancelling")]

    # -- every second ----------------------------------------------------------
    def poll(self):
        live = set()
        for r in self.jobs():
            ident, p = r["id"], r["profile"]
            job = (r.get("remote") or {}).get("job_id")
            if p.get("backend") == "slurm" and job:
                # Queued, or ended with results still downloading: nothing to read.
                if (r.get("remote") or {}).get("slurm_state") not in (None, "RUNNING") or ident in self.finished:
                    continue
                live.add(ident)
                if ident not in self.nodes:
                    self.ask_node(ident, p["host"], job)
                elif ident not in self.streams and time.time() - self._opened.get(ident, 0) >= NODE_RETRY_S:
                    self.stream(ident, p["host"], job)          # a dropped stream reconnects at most every 15 s
            elif r.get("pid"):
                live.add(ident)
                self.local(r)
        for ident in list(self.streams):                # jobs that ended
            if ident not in live:
                self.stop(ident)
        for ident in list(self.latest):
            if ident not in live:
                self.latest.pop(ident, None)
        self.updated.emit()

    # -- cluster -----------------------------------------------------------------
    def ask_node(self, ident, host, job):
        now = time.time()
        if now - self._asked.get(ident, 0) < NODE_RETRY_S:
            return
        self._asked[ident] = now
        proc = QtCore.QProcess(self)
        proc.finished.connect(lambda *_a: self.node_answer(ident, proc))
        proc.start("ssh", ["-o", "ConnectTimeout=5", "-o", "BatchMode=yes", host,
                           "sacct -n -P -X -j %s -o NodeList,AllocCPUS,ReqMem,State" % job])

    def node_answer(self, ident, proc):
        out = bytes(proc.readAllStandardOutput()).decode("utf8", "replace").split()
        proc.deleteLater()
        f = out[0].split("|") if out else []
        if len(f) >= 4 and f[3].startswith("RUNNING") and f[0] and f[0] != "None assigned":
            self.nodes[ident] = (f[0], int(f[1]) if f[1].isdigit() else None, _gb(f[2]))
        elif len(f) >= 4 and f[3] and not f[3].startswith(("PENDING", "REQUEUED", "RESIZING", "SUSPENDED")):
            self.finished.add(ident)                    # ended: never reconnect

    def stream(self, ident, host, job):
        node = self.nodes[ident][0]
        proc = QtCore.QProcess(self)
        proc.setProcessChannelMode(QtCore.QProcess.SeparateChannels)
        proc.readyReadStandardOutput.connect(lambda: self.read(ident, proc))
        proc.finished.connect(lambda *_a: self.ended(ident, proc))
        self.streams[ident] = proc
        self._opened[ident] = time.time()
        proc.start("ssh", ["-o", "ConnectTimeout=5", "-o", "BatchMode=yes", "-o", "ServerAliveInterval=15", host,
                           "ssh -o BatchMode=yes -o ConnectTimeout=5 %s bash -s" % node])
        proc.write((NODE_LOOP % dict(job=job)).encode())
        proc.closeWriteChannel()

    def read(self, ident, proc):
        node, cores, alloc = self.nodes.get(ident, ("", None, None))
        while proc.canReadLine():
            f = bytes(proc.readLine()).decode("utf8", "replace").split()
            if f == ["END"]:
                self.finished.add(ident)
            if len(f) != 8 or f[0] != "S":
                continue
            try:
                t, usage, inactive, peak, limit, cpu_ns, ncpu = (float(x or 0) for x in f[1:])
            except ValueError:
                continue
            d = self.latest.get(ident) or {}
            busy = None
            if d.get("_cpu_ns") is not None and t - d["_t"] > 0.2 and cpu_ns >= d["_cpu_ns"]:
                busy = (cpu_ns - d["_cpu_ns"])/1e9/(t - d["_t"])
            self.latest[ident] = dict(node=node, cores=int(ncpu) or cores, alloc_gb=limit/1e9 or alloc,
                                      rss_gb=(usage - inactive)/1e9, peak_gb=peak/1e9, cores_busy=busy,
                                      t=time.time(), _t=t, _cpu_ns=cpu_ns)
        self.updated.emit()

    def ended(self, ident, proc):
        if self.streams.get(ident) is proc:
            self.streams.pop(ident)
            self.nodes.pop(ident, None)                 # Slurm closes it when the job ends: ask sacct again
        try:
            proc.deleteLater()
        except RuntimeError:                            # already gone (monitor torn down)
            pass

    def stop(self, ident):
        proc = self.streams.pop(ident, None)
        if proc is not None:
            proc.kill()
            proc.waitForFinished(1000)
        self.nodes.pop(ident, None)

    def shutdown(self):
        """Stop polling and close every stream (before a reload replaces this monitor)."""
        self.timer.stop()
        for ident in list(self.streams):
            self.stop(ident)

    # -- local -------------------------------------------------------------------
    def local(self, r):
        pid = int(r["pid"])
        try:
            rss = cpu = 0.0
            for p in self._tree(pid):
                rss += self._pss_gb(p)
                with open("/proc/%d/stat" % p) as f:
                    fields = f.read().rsplit(")", 1)[1].split()
                cpu += (int(fields[11]) + int(fields[12]))/CLK
        except (OSError, ValueError, IndexError):
            self.latest.pop(r["id"], None)
            return
        now = time.time()
        prev = self._cpu.get(r["id"])
        busy = max(0.0, (cpu - prev[0])/(now - prev[1])) if prev and now - prev[1] > 0.2 else None
        self._cpu[r["id"]] = (cpu, now)
        self.latest[r["id"]] = dict(node="this machine", cores=os.cpu_count(), alloc_gb=None, rss_gb=rss,
                                    peak_gb=max(rss, (self.latest.get(r["id"]) or {}).get("peak_gb") or 0),
                                    cores_busy=busy, t=now)

    @staticmethod
    def _pss_gb(pid):
        try:
            with open("/proc/%d/smaps_rollup" % pid) as f:
                m = re.search(r"^Pss:\s+(\d+)", f.read(), re.M)
            if m:
                return int(m.group(1))*1e-6
        except OSError:
            pass
        with open("/proc/%d/status" % pid) as f:
            m = re.search(r"VmRSS:\s+(\d+)", f.read())
        return int(m.group(1))*1e-6 if m else 0.0

    @staticmethod
    def _tree(pid):
        kids = {}
        for name in os.listdir("/proc"):
            if name.isdigit():
                try:
                    with open("/proc/%s/stat" % name) as f:
                        ppid = int(f.read().rsplit(")", 1)[1].split()[1])
                    kids.setdefault(ppid, []).append(int(name))
                except (OSError, ValueError, IndexError):
                    pass
        out, todo = [], [pid]
        while todo:
            p = todo.pop()
            out.append(p)
            todo += kids.get(p, [])
        return out

    # -- reading -----------------------------------------------------------------
    def sample(self, ident):
        d = self.latest.get(ident)
        if d is None:
            return None
        return dict({k: v for k, v in d.items() if not k.startswith("_")}, age_s=time.time() - d["t"])

    def line(self, ident):
        d = self.sample(ident)
        if d is None:
            return ""
        parts = []
        if d.get("rss_gb") is not None:
            ram = ("RAM %.0f" if d["rss_gb"] >= 10 else "RAM %.1f") % d["rss_gb"]
            ram += (" / %.0f GB" % d["alloc_gb"]) if d.get("alloc_gb") else " GB"
            if d.get("peak_gb") and d["peak_gb"] > 1.05*d["rss_gb"]:
                ram += " (peak %.0f)" % d["peak_gb"]
            parts.append(ram)
        if d.get("cores_busy") is not None:
            parts.append("%.1f of %s cores busy" % (d["cores_busy"], d.get("cores") or "?"))
        if d.get("node"):
            parts.append(d["node"])
        if d["age_s"] > 5:
            parts.append("%d s ago" % d["age_s"])
        return " · ".join(parts)


_MONITOR = None


def monitor():
    global _MONITOR
    if _MONITOR is None:
        _MONITOR = Monitor()
    return _MONITOR
