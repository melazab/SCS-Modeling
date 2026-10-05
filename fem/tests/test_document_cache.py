from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from scs_document_cache import same_surface, same_params


class DocumentCacheTests(unittest.TestCase):
    def setUp(self):
        self.v = np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]])
        self.t = np.array([[0,1,2],[0,2,3]])

    def test_vertex_and_triangle_numbering(self):
        order = np.array([2,0,3,1]); inverse = np.argsort(order)
        self.assertTrue(same_surface((self.v[order],inverse[self.t[::-1]]), (self.v,self.t)))

    def test_planar_diagonal_and_stl_rounding(self):
        alternate = np.array([[0,1,3],[1,2,3]])
        self.assertTrue(same_surface((self.v+1e-6,alternate), (self.v,self.t)))

    def test_warped_patch_and_real_movement_rejected(self):
        warped = self.v.copy(); warped[2,2] = .01
        alternate = np.array([[0,1,3],[1,2,3]])
        self.assertFalse(same_surface((warped,alternate),(warped,self.t)))
        self.assertFalse(same_surface((self.v+.001,self.t),(self.v,self.t)))

    def test_missing_face_rejected(self):
        self.assertFalse(same_surface((self.v,self.t[:1]),(self.v,self.t)))

    def test_float_ui_rounding_not_resolution_changes(self):
        self.assertTrue(same_params({'FIELDS':[('bone',.7499999999998)]},
                                    {'FIELDS':[['bone',.75]]}))
        self.assertFalse(same_params({'H_MAX':.8},{'H_MAX':.81}))


if __name__ == '__main__':
    unittest.main()
