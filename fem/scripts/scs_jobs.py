"""Shared execution profiles and persistent job records, independent of FreeCAD UI.

A remote profile is an SSH host plus sbatch options; it can be checked, but
remote submission from the panel is not implemented (fem/hpc/ scripts run jobs).
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
    return {'local': local, 'case-hpc': dict(name='Case HPC', backend='slurm', host='case-hpc',
            sbatch_options='-n 16 --mem=96gb --time=04:00:00')}


# The repository mirror on the cluster, relative to the remote $HOME, and the
# worker environment there (Elmer, Python venv); see fem/hpc/scs-env.sh.
REMOTE_ROOT = 'scs/SCS-Modeling'
REMOTE_ENV = '$HOME/%s/fem/hpc/scs-env.sh' % REMOTE_ROOT
SAFE_OPTION = re.compile(r'[A-Za-z0-9_.,:=@%/+-]+')


def sbatch_options(profile):
    """The profile's sbatch options as a token list. Older profiles stored
    separate account/partition/CPU/memory/time fields; they convert here."""
    if 'sbatch_options' in profile:
        text = profile['sbatch_options']
    else:
        text = '-n %d --mem=%dgb --time=%s' % (int(profile.get('cpus', 16)),
                                              int(profile.get('memory_gb', 64)),
                                              profile.get('time_limit', '04:00:00'))
        if profile.get('account'):
            text += ' -A ' + profile['account']
        if profile.get('partition'):
            text += ' -p ' + profile['partition']
    try:
        tokens = shlex.split(text)
    except ValueError as exc:
        raise ValueError('sbatch options: %s' % exc)
    if not tokens or not tokens[0].startswith('-') or not all(SAFE_OPTION.fullmatch(t) for t in tokens):
        raise ValueError('sbatch options must be plain options such as -n 24 --mem=140gb '
                         '--time=06:00:00 (no quotes, spaces in values or shell characters).')
    return tokens


def validate_profile(profile):
    p = copy.deepcopy(profile)
    if p['backend'] == 'local':
        p['mesh_threads'] = max(1, min(os.cpu_count() or 1, int(p['mesh_threads'])))
        p['solve_cpus'] = max(1, min(physical_core_limit(), int(p['solve_cpus'])))
    elif p['backend'] == 'slurm':
        if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.@-]*', p['host']):
            raise ValueError('Use an SSH host alias such as case-hpc; no shell commands.')
        p = dict(name=p['name'], backend='slurm', host=p['host'],
                 sbatch_options=' '.join(sbatch_options(p)))
        return p
    else:
        raise ValueError('Unknown execution backend.')
    p['memory_gb'] = float(p['memory_gb'])
    if not 0 <= p['memory_gb'] < 1e6:
        raise ValueError('Memory must be a finite nonnegative GB value.')
    return p


def requested(tokens):
    """(tasks x cpus-per-task, memory MB or None, mem-per-cpu MB or None, account, partition)."""
    value = {}
    flags = {'-n': 'ntasks', '--ntasks': 'ntasks', '-c': 'cpus', '--cpus-per-task': 'cpus',
             '--mem': 'mem', '--mem-per-cpu': 'mem_per_cpu', '-A': 'account', '--account': 'account',
             '-p': 'partition', '--partition': 'partition'}
    i = 0
    while i < len(tokens):
        name, eq, rest = tokens[i].partition('=')
        if name in flags:
            if not eq and i + 1 < len(tokens):
                i += 1
                rest = tokens[i]
            value[flags[name]] = rest
        i += 1

    def mb(text):
        m = re.fullmatch(r'(\d+)([KMGT]?)B?', (text or '').upper())
        return None if not m else int(m.group(1)) * {'K': 1 / 1024, '': 1, 'M': 1, 'G': 1024, 'T': 1024 ** 2}[m.group(2)]
    cpus = int(value.get('ntasks', 1)) * int(value.get('cpus', 1))
    return cpus, mb(value.get('mem')), mb(value.get('mem_per_cpu')), value.get('account', ''), value.get('partition', '')


def allocation_gb(profile):
    """Memory the profile's sbatch options request, in GB, or None if unset."""
    cpus, mem, mem_per_cpu, _, _ = requested(sbatch_options(profile))
    mb = mem if mem else (mem_per_cpu * cpus if mem_per_cpu else None)
    return None if mb is None else mb / 1024


def ssh_probe_args(profile):
    """Read-only check of the cluster as the jobs will see it.

    Loads the worker environment and reports each tool as OK/MISSING. It then
    reads the partition's MaxMemPerCPU and the account's group CPU limit:
    SLURM silently raises the CPU count to cover --mem, and a job over the
    group limit pends forever, which sbatch --test-only does not reveal.
    Nothing is submitted. Exits 0 only when everything needed is present."""
    p = validate_profile(profile)
    tokens = sbatch_options(p)
    _, _, _, account, partition = requested(tokens)
    command = (
        'missing=0; echo "OK host $(hostname)"; '
        'if [ -f %(env)s ]; then source %(env)s >/dev/null 2>&1; echo "OK worker environment %(env)s"; '
        'else echo "MISSING worker environment %(env)s"; missing=1; fi; '
        'for c in sbatch squeue scancel sacct ElmerGrid ElmerSolver_mpi mpirun; do '
        'p=$(command -v $c) && echo "OK $c $p" || { echo "MISSING $c"; missing=1; }; done; '
        'if [ -n "$SCS_PY" ] && "$SCS_PY" -c "import numpy, scipy, gmsh, vtk" 2>/dev/null; '
        'then echo "OK python $SCS_PY (numpy, scipy, gmsh, vtk)"; '
        'else echo "MISSING python environment ${SCS_PY:-unset}"; missing=1; fi; '
        'part=%(part)s; [ -n "$part" ] || part=$(scontrol show partition -o | sed -n "s/^PartitionName=\\([^ ]*\\).*Default=YES.*/\\1/p"); '
        'acct=%(acct)s; [ -n "$acct" ] || acct=$(sacctmgr -nP show user $USER format=defaultaccount); '
        'echo "LIMIT partition=$part account=$acct '
        'maxmempercpu=$(scontrol show partition "$part" -o | grep -o "MaxMemPerCPU=[0-9]*" | cut -d= -f2) '
        'groupcpu=$(sacctmgr -nP show assoc account="$acct" user= format=grptres | grep -o "cpu=[0-9]*" | cut -d= -f2) '
        'groupused=$(squeue -h -A "$acct" -t RUNNING -o %%C | paste -sd+ | bc 2>/dev/null) '
        'cpus=$(sinfo -h -p "$part" -o %%C)"; '
        'exit $missing' % dict(env=REMOTE_ENV, part=shlex.quote(partition), acct=shlex.quote(account)))
    # Fixed commands; the only inputs are a validated alias and option tokens.
    return ['-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', '-o',
            'StrictHostKeyChecking=yes', p['host'], 'bash -lc ' + shlex.quote(command)]


def describe_probe(code, output, profile):
    """(one-line verdict, details) for ssh_probe_args() output."""
    lines = output.strip().splitlines()
    missing = [l[len('MISSING '):] for l in lines if l.startswith('MISSING ')]
    if code == 255:
        return 'Could not connect (exit 255). Is the VPN up and the host key trusted?', output
    if missing:
        return 'Connected, but missing: ' + ', '.join(missing), output
    if code != 0:
        return 'Check failed with exit code %d.' % code, output
    limit = dict(kv.split('=', 1) for l in lines if l.startswith('LIMIT ')
                 for kv in l[len('LIMIT '):].split() if '=' in kv)
    cpus, mem, mem_per_cpu, _, _ = requested(sbatch_options(profile))
    per_cpu = int(limit['maxmempercpu']) if limit.get('maxmempercpu', '').isdigit() else None
    group = int(limit['groupcpu']) if limit.get('groupcpu', '').isdigit() else None
    # SLURM raises CPUs per task until memory per CPU fits, so the total
    # grows in whole multiples: -n 24 --mem=200gb became 48 CPUs, not 34.
    per_task = mem_per_cpu or (mem / cpus if mem else None)
    granted = cpus * max(1, -(-int(per_task) // per_cpu)) if per_cpu and per_task else cpus
    where = 'account %s, partition %s' % (limit.get('account', '?'), limit.get('partition', '?'))
    if group and granted > group:
        verdict = ('These options would never start: they need %d CPUs and %s allows %d in total.'
                   % (granted, where, group))
        if granted > cpus:
            verdict += ' The memory request is what raises it (at most %d MB per CPU); ask for <= %d GB.' % (
                per_cpu, cpus * per_cpu // 1024)
    elif granted > cpus:
        verdict = ('Ready, but SLURM will charge %d CPUs, not %d: the memory request exceeds %d MB per CPU (%s).'
                   % (granted, cpus, per_cpu, where))
    else:
        verdict = 'Ready: %d CPUs on %s%s.' % (granted, where, ', group limit %d' % group if group else '')
    # Live load: the group's running CPUs, and sinfo's allocated/idle/other/total.
    used = int(limit['groupused']) if limit.get('groupused', '').isdigit() else 0
    if group and used and granted <= group < used + granted:
        verdict += ' Your group is already running %d CPUs, so it will wait for those to finish.' % used
    load = limit.get('cpus', '').split('/')
    if len(load) == 4:
        verdict += ' %s now has %s of %s CPUs idle.' % (limit.get('partition', 'The partition'), load[1], load[3])
    return verdict, output.strip()


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

    def clear(self, ident):
        # Archive, never delete: the record and its copied log stay under cleared/.
        cleared = self.root / 'cleared'
        cleared.mkdir(parents=True, exist_ok=True)
        for name in (ident+'.json', ident+'.log'):
            if (self.root / name).exists():
                (self.root / name).replace(cleared / name)

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
