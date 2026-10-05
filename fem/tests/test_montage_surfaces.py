"""Coverage of anatomy and disconnected bodies in potential surfaces."""
import sys
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import scs_montage as M

class SurfaceTests(unittest.TestCase):
    def test_all_classified_anatomy_is_visualizable(self):
        paths = M.surface_paths()
        self.assertEqual(set(paths), set(M.C.TISSUE_BODIES))
        for tissue, bodies in M.C.TISSUE_BODIES.items():
            self.assertEqual(tuple(bodies), paths[tissue])
            self.assertIn(tissue, M.TISSUE_LABELS)

    def test_disconnected_bodies_get_correct_connectivity_and_values(self):
        vertices = np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        triangles = np.array([[0,1,2]])
        with patch.object(M.stlio, 'read_stl', side_effect=[(vertices,triangles),(vertices+3,triangles)]), \
             patch.object(M, 'montage_at_points', side_effect=lambda currents, points, **kw: points[:,0]):
            v, t, values = M.sample_tissue('root', [1,-1], surface_path=('a.stl','b.stl'))
        np.testing.assert_array_equal(t, [[0,1,2],[3,4,5]])
        np.testing.assert_array_equal(values, v[:,0])
        np.testing.assert_array_equal(v[t[1]], vertices+3)

if __name__ == '__main__':
    unittest.main()
