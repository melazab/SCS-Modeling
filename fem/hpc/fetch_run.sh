#!/bin/bash
# Copy a finished Pioneer run back to this workstation:
#   fem/hpc/fetch_run.sh <run_name> [local_name]
# <run_name> is the directory under ~/scs/SCS-Modeling/fem/out/lead_runs/ on
# Pioneer; it lands in fem/out/lead_runs/<local_name> here. Only the published
# artifacts come back (mesh.npz, preview.vtp, solution*.npz, their manifests,
# reports, params and the lead STLs). Elmer's scratch (elmer*/, tens of GB),
# mesh.msh and the size field stay on the cluster.
set -euo pipefail
NAME=$1; LOCAL=${2:-$1}
HERE=$(cd "$(dirname "$0")/../.." && pwd)
SRC="scs/SCS-Modeling/fem/out/lead_runs/$NAME"
DEST="$HERE/fem/out/lead_runs/$LOCAL"
mkdir -p "$DEST"
rsync -a --partial --info=progress2 -e "ssh -o ConnectTimeout=5" \
  --include='mesh.npz*' --include='preview.vtp*' --include='solution*.npz*' \
  --include='*_report.json' --include='params.json' --include='lead/***' --exclude='*' \
  "case-hpc:$SRC/" "$DEST/"
# Paths inside the reports point at the cluster; repoint them here. They are
# not part of any provenance signature.
python3 - "$DEST" "/home/[^/]+/$SRC" <<'PY'
import json, re, sys
from pathlib import Path
dest, remote = Path(sys.argv[1]), re.compile(sys.argv[2])
for name in ('params.json', 'mesh_report.json', 'solve_report.json'):
    path = dest / name
    if path.exists():
        path.write_text(remote.sub(str(dest), path.read_text()))
PY
# Re-hash what arrived against its manifest; a truncated copy must not pass.
python3 - "$DEST" <<'PY'
import hashlib, json, sys
from pathlib import Path
for manifest in sorted(Path(sys.argv[1]).glob('*.manifest.json')):
    path = manifest.with_name(manifest.name[:-len('.manifest.json')])
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    ok = h.hexdigest() == json.loads(manifest.read_text())['sha256']
    print('%-22s %s' % (path.name, 'sha256 OK' if ok else 'SHA256 MISMATCH'))
    if not ok:
        sys.exit(1)
PY
echo "fetched into $DEST"
