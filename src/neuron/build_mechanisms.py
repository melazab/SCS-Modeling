"""Compile the pinned, unmodified channel sources into a private build directory."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BUILD = HERE / 'build'


def build():
    import neuron
    records = json.loads((HERE / 'vendor/sources.json').read_text())
    for entry in records:
        path = HERE / 'vendor' / entry['path']
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            raise RuntimeError(f'Vendored source changed: {path}')
    sources = sorted((HERE / 'vendor').rglob('*.mod'))
    signature = hashlib.sha256((''.join(p.read_text() for p in sources)
                                + neuron.__version__ + sys.version).encode()).hexdigest()
    stamp = BUILD / 'signature.txt'
    libraries = list(BUILD.glob('*/libnrnmech.so'))
    if stamp.exists() and stamp.read_text() == signature and libraries:
        return libraries[0]
    BUILD.mkdir(exist_ok=True)
    mods = BUILD / 'mods'
    mods.mkdir(exist_ok=True)
    for source in sources:
        shutil.copyfile(source, mods / source.name)
    exe = Path(sys.executable).parent / 'nrnivmodl'
    result = subprocess.run([str(exe), str(mods)], cwd=BUILD, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (BUILD / 'compile.log').write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(result.stdout[-6000:])
    stamp.write_text(signature)
    return next(BUILD.glob('*/libnrnmech.so'))


if __name__ == '__main__':
    print(build())
