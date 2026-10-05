"""Physical cores available to the local solve, respecting Linux affinity."""
import os
from pathlib import Path


def physical_core_limit():
    allowed = os.sched_getaffinity(0) if hasattr(os, 'sched_getaffinity') else range(os.cpu_count() or 1)
    try:
        cores = set()
        for cpu in allowed:
            topology = Path('/sys/devices/system/cpu/cpu%d/topology' % cpu)
            cores.add((topology.joinpath('physical_package_id').read_text().strip(),
                       topology.joinpath('core_id').read_text().strip()))
        return max(1, len(cores))
    except OSError:
        try:
            import psutil
            return max(1, min(psutil.cpu_count(logical=False) or 1, len(allowed)))
        except ImportError:
            return 1  # unknown topology: do not promise unsupported MPI slots
