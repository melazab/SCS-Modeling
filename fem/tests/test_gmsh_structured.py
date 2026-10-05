"""Exercise the real Structured field and HXT in an isolated native process."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class StructuredMeshingTests(unittest.TestCase):
    def test_parallel_volume_with_serial_boundary(self):
        scripts = Path(__file__).resolve().parents[1] / 'scripts'
        code = r'''
import sys
import numpy as np
import gmsh
from build_mesh import configure_gmsh_threads
path = sys.argv[1]
with open(path, 'w') as f:
    f.write('0 0 0\n0.1 0.1 0.1\n21 21 21\n')
    np.full(21**3, .25).tofile(f, sep='\n')
    f.write('\n')
gmsh.initialize()
try:
    configure_gmsh_threads(gmsh, 20)
    gmsh.model.occ.addBox(0, 0, 0, 2, 2, 2)
    gmsh.model.occ.synchronize()
    field = gmsh.model.mesh.field.add('Structured')
    gmsh.model.mesh.field.setString(field, 'FileName', path)
    gmsh.model.mesh.field.setNumber(field, 'TextFormat', 1)
    gmsh.model.mesh.field.setAsBackgroundMesh(field)
    gmsh.option.setNumber('Mesh.Algorithm3D', 10)
    for name in ('MeshSizeFromPoints', 'MeshSizeFromCurvature', 'MeshSizeExtendFromBoundary'):
        gmsh.option.setNumber('Mesh.'+name, 0)
    gmsh.model.mesh.generate(3)
    assert gmsh.option.getNumber('Mesh.MaxNumThreads3D') == 20
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    tet_tags, conn = gmsh.model.mesh.getElementsByType(4)
    assert len(tet_tags) > 100
    order = np.argsort(tags)
    p = xyz.reshape(-1,3)[order][np.searchsorted(tags[order],conn)].reshape(-1,4,3)
    volumes = np.linalg.det(p[:,1:] - p[:,:1]) / 6
    assert np.all(volumes > 0), volumes.min()
    assert abs(volumes.sum() - 8) < 1e-8, volumes.sum()
finally:
    gmsh.finalize()
'''
        with tempfile.TemporaryDirectory() as tmp:
            import os
            env = dict(os.environ, PYTHONPATH=str(scripts), OMP_NUM_THREADS='20')
            result = subprocess.run([sys.executable, '-c', code, str(Path(tmp)/'field.dat')],
                                    env=env, capture_output=True, text=True, timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
