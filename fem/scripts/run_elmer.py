"""Batch-workflow entry point for the same Elmer backend as the visualizer.

The live-lead GUI uses solve_lead.py with run-specific provenance. This entry
point serves run_all.sh's configured, frozen lead and writes fem/out results.
"""
import argparse
import os
import numpy as np
import artifacts
import config as C
from elmer_backend import solve
from solve_resources import Resources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpus', type=int, default=int(os.environ.get('OMP_NUM_THREADS', '1')))
    parser.add_argument('--memory-gb', type=float, default=float(os.environ.get('SCS_SOLVE_MEMORY_GB', '0')))
    args = parser.parse_args()
    if hasattr(os, 'setsid') and os.getpgrp() != os.getpid():
        os.setsid()
    background = os.environ.get('SCS_BACKGROUND', '0') == '1'
    with Resources(args.cpus, args.memory_gb) as resources:
        with np.load(os.path.join(C.OUT, 'mesh.npz')) as mesh:
            nodes, tets = mesh['nodes'], mesh['tets']
        result = solve(nodes, tets, C.CONTACT_STL, C.INSULATOR_STL, background,
                       os.path.join(C.OUT, 'elmer_batch_bg' if background else 'elmer_batch'),
                       resources, lambda value: print('PROGRESS %d' % value, flush=True))
        artifacts.atomic_npz(artifacts.result_path('solution.npz'), **result)


if __name__ == '__main__':
    main()
