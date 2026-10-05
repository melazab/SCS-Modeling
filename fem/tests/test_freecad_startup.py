"""Exercise the addon using FreeCAD's exec-based loader semantics."""
import os
from pathlib import Path
import subprocess
import sys
import unittest


class StartupTests(unittest.TestCase):
    def test_loader_file_does_not_determine_repo(self):
        addon = Path(__file__).resolve().parents[1] / 'freecad_startup'
        code = r'''
import pathlib, sys, types
addon = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(addon))
# Isolate this bootstrap test from FreeCAD and numerical dependencies.
sys.modules['scs_montage_feature'] = types.ModuleType('scs_montage_feature')
for filename in (None, '/unrelated/FreeCADInit.py'):
    expected = str(addon.parent / 'scripts')
    sys.path[:] = [p for p in sys.path if p != expected]
    namespace = {} if filename is None else {'__file__': filename}
    exec(compile((addon / 'Init.py').read_text(), str(addon / 'Init.py'), 'exec'), namespace)
    assert sys.path[0] == expected, sys.path
'''
        env = dict(os.environ)
        env.pop('SCS_MODELING_REPO', None)
        subprocess.run([sys.executable, '-c', code, str(addon)], env=env, check=True)


if __name__ == '__main__':
    unittest.main()
