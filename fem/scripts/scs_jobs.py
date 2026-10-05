"""Shared execution profiles and persistent job records, independent of FreeCAD UI.

Remote profiles can be configured/probed; remote submission is not implemented.
Never silently fall back to local execution when a remote profile is selected.
"""
import copy
import json
import os
from pathlib import Path
import re
import shlex
import time
import uuid

from solve_cpu_limit import physical_core_limit

ROOT = Path(__file__).resolve().parents[1] / 'out' / 'jobs'
TERMINAL = {'completed', 'failed', 'cancelled', 'unverified'}
TRANSITIONS = {
    'starting': {'running', 'failed', 'cancelling', 'completed'},
    'running': {'completed', 'failed', 'cancelling'},
    'cancelling': {'cancelled', 'failed', 'completed'},
}


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(data, indent=2, sort_keys=True))
    temp.replace(path)


def defaults():
    local = dict(name='Local workstation', backend='local',
                 mesh_threads=os.cpu_count() or 1,
                 solve_cpus=min(4, physical_core_limit()), memory_gb=0.)
    try:
        legacy = json.loads((ROOT.parent / 'solve_settings.json').read_text())
        local.update(solve_cpus=min(physical_core_limit(), max(1,int(legacy['cpus']))),
                     memory_gb=float(legacy.get('memory_gb', 0.)))
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return {'local': local, 'case-hpc': dict(name='Case HPC', backend='slurm',
            host='case-hpc', account='', partition='', directory='',
            cpus=16, memory_gb=64., time_limit='04:00:00')}


def validate_profile(profile):
    p = copy.deepcopy(profile)
    if p['backend'] == 'local':
        p['mesh_threads'] = max(1, min(os.cpu_count() or 1, int(p['mesh_threads'])))
        p['solve_cpus'] = max(1, min(physical_core_limit(), int(p['solve_cpus'])))
    elif p['backend'] == 'slurm':
        if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.@-]*', p['host']):
            raise ValueError('Use an SSH host alias such as case-hpc; no shell commands.')
        if not re.fullmatch(r'(?:\d+-)?\d{1,3}:\d{2}:\d{2}', p['time_limit']):
            raise ValueError('Time limit must be HH:MM:SS or D-HH:MM:SS.')
        p['cpus'] = max(1, int(p['cpus']))
        for key in ('account','partition','directory'):
            if any(c in p[key] for c in '\r\n\x00'):
                raise ValueError('Invalid %s.' % key)
    else:
        raise ValueError('Unknown execution backend.')
    p['memory_gb'] = float(p['memory_gb'])
    if not 0 <= p['memory_gb'] < 1e6:
        raise ValueError('Memory must be a finite nonnegative GB value.')
    return p


def ssh_probe_args(profile):
    p = validate_profile(profile)
    command = ('hostname; command -v sbatch; command -v squeue; '
               'command -v scancel; command -v sacct; '
               'command -v python3; command -v ElmerSolver_mpi; '
               'command -v ElmerGrid')
    # Fixed read-only commands. Login shell loads the cluster module environment.
    return ['-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', '-o',
            'StrictHostKeyChecking=yes', p['host'], 'bash -lc ' + shlex.quote(command)]


class Store:
    def __init__(self, root=ROOT):
        self.root = Path(root)
        self.profiles = defaults()
        if (self.root / 'profiles.json').exists():
            self.profiles.update(json.loads((self.root / 'profiles.json').read_text()))

    def save_profile(self, key, profile):
        self.profiles[key] = validate_profile(profile)
        atomic_json(self.root / 'profiles.json', self.profiles)

    def new_job(self, document, task, profile, args, environment, log_path):
        record = dict(id=uuid.uuid4().hex, document=document, task=task,
                      profile=copy.deepcopy(profile), args=list(args),
                      environment=dict(environment or {}), state='starting',
                      progress=0, created=time.time(), updated=time.time(), log_path=log_path)
        self.write(record)
        return record

    def write(self, record):
        atomic_json(self.root / (record['id']+'.json'), record)

    def transition(self, record, state, **fields):
        old = record['state']
        if state != old and state not in TRANSITIONS.get(old, set()):
            raise ValueError('Invalid job transition %s -> %s' % (old, state))
        record.update(fields, state=state, updated=time.time())
        self.write(record)

    def history(self):
        records=[]
        for path in self.root.glob('*.json'):
            if path.name=='profiles.json': continue
            try:
                r=json.loads(path.read_text())
                if 'id' in r: records.append(r)
            except (OSError, ValueError):
                continue
        return sorted(records,key=lambda r:r['created'],reverse=True)
