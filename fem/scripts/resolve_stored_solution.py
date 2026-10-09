"""Re-solve a published solution's exact classified mesh and compare fields.

CLI: resolve_stored_solution.py <solution_bg.npz> <work_dir> [--cpus N]
     [--memory-gb G] [--contacts 1 2 ...]

Feeds the stored nodes, tets, labels, conductivities and contact order straight
to elmer_backend.solve_prepared(), so the only thing that differs from the run
that produced the file is the Elmer build/MPI layout. Used to validate the
Pioneer build against the workstation. Classification is deliberately skipped:
a config.py change must not be mistaken for a solver difference.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elmer_backend import ELMER_BIN, solve_prepared
from solve_resources import Resources


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('solution')
    p.add_argument('work_dir')
    p.add_argument('--cpus', type=int, default=1)
    p.add_argument('--memory-gb', type=float, default=0.)
    p.add_argument('--contacts', type=int, nargs='*', help='subset to solve; default all stored contacts')
    args = p.parse_args()

    ref = np.load(args.solution)
    order = [str(x) for x in ref['order']]
    stored = [int(c) for c in ref['contact_ids']]
    contacts = args.contacts or stored
    missing = sorted(set(contacts) - set(stored))
    if missing:
        raise SystemExit('contacts %s are not in the stored solution' % missing)
    print('Elmer: %s' % ELMER_BIN, flush=True)
    t0 = time.time()
    with Resources(args.cpus, args.memory_gb) as resources:
        nodes, tets = ref['nodes'], ref['tets']
        result = solve_prepared(nodes, tets, ref['label'], ref['sigma'], order, contacts,
                                args.work_dir, resources,
                                lambda pct: print('PROGRESS %d' % pct, flush=True))
    if int(result['pinned_node']) != int(ref['pinned_node']):
        raise SystemExit('pinned node differs (%d vs %d); fields are not comparable'
                         % (result['pinned_node'], ref['pinned_node']))

    phi_ref = ref['phi']
    rows = []
    for pos, contact in enumerate(contacts):
        a, b = result['phi'][pos], phi_ref[stored.index(contact)]
        span = float(b.max() - b.min())
        diff = np.abs(a - b)
        rows.append(dict(contact=contact, span_v=span, max_abs_v=float(diff.max()),
                         max_rel_of_span=float(diff.max() / span),
                         rms_rel_of_span=float(np.sqrt(np.mean(diff ** 2)) / span),
                         residual=float(result['relative_residuals'][pos])))
        print('contact %2d: span %.4g V, max |diff| %.3e V (%.2e of span), rms %.2e of span'
              % (contact, span, diff.max(), diff.max() / span, rows[-1]['rms_rel_of_span']), flush=True)
    report = dict(solution=os.path.abspath(args.solution), elmer_bin=str(ELMER_BIN),
                  cpus=args.cpus, reference_cpus=int(ref['cpus']), nodes=len(nodes), tets=len(tets),
                  elapsed_s=round(time.time() - t0, 1), peak_rss_gb=resources.peak_gb,
                  worst_rel_of_span=max(r['max_rel_of_span'] for r in rows), contacts=rows)
    with open(os.path.join(args.work_dir, 'comparison.json'), 'w') as fh:
        json.dump(report, fh, indent=2)
    print('COMPARE_JSON ' + json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
