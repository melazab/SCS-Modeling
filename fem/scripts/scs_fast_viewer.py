"""Self-contained WebGL potential viewer: fixed scale, real timing slowed uniformly."""
import json
from pathlib import Path
import numpy as np
import scs_waveform as W


def export(path, surfaces, metadata):
    W.validate(metadata['frequency_hz'], metadata['width_us'], metadata['gap_us'])
    W.require_balanced(metadata['currents_mA'])
    clamp = float(metadata['clamp_v'])
    if not np.isfinite(clamp) or clamp <= 0:
        raise ValueError('Positive fixed color clamp required')
    data = dict(metadata=metadata, surfaces=[])
    for name, (vertices, triangles, values) in surfaces.items():
        v, t, a = np.asarray(vertices), np.asarray(triangles), np.asarray(values)
        if not np.isfinite(v).all() or not np.isfinite(a).all():
            raise ValueError('Cannot export nonfinite surface values: '+name)
        data['surfaces'].append(dict(name=name, xyz=v[t].reshape(-1, 3).tolist(),
                                     values=a[t].ravel().tolist()))
    payload = json.dumps(data, allow_nan=False).replace('<', '\\u003c')
    html = Path(__file__).with_name('scs_fast_viewer.html').read_text()
    Path(path).write_text(html.replace('/*SCS_DATA*/', 'const data = '+payload+';'))
    Path(path).with_suffix('.json').write_text(json.dumps(metadata, indent=2))
    return str(path)
