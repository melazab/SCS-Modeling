from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from mesh_geometry import Extent


class CropWindowTests(unittest.TestCase):
    def setUp(self):
        # A canal that drifts 0.2 mm in y per mm of z, 0..200 mm, and a lead at z 90..130.
        z = np.linspace(0, 200, 2001)
        canal = np.column_stack([np.full_like(z, 50.), 10 + 0.2 * z, z])
        lead = np.column_stack([np.full(41, 50.), 10 + 0.2 * np.linspace(90, 130, 41), np.linspace(90, 130, 41)])
        self.clouds = {'canal': canal, 'lead': lead}
        pts = np.vstack([canal, lead])
        self.extent = Extent(self.clouds, np.array([pts.min(0), pts.max(0)]))

    def test_zero_crop_is_the_whole_model(self):
        lo, hi = self.extent.box(0, 1.0)
        np.testing.assert_allclose(lo, [49, 9, -1])
        np.testing.assert_allclose(hi, [51, 51, 201])

    def test_crop_cuts_z_and_hugs_the_anatomy_inside(self):
        lo, hi = self.extent.box(20, 1.0)
        self.assertAlmostEqual(lo[2], 70)          # 20 mm below the lead, no margin on the cut
        self.assertAlmostEqual(hi[2], 150)
        # y covers only the canal between z 70 and 150 (24..40), plus the margin and one slab
        self.assertGreater(lo[1], 22)
        self.assertLess(hi[1], 42)

    def test_crop_never_extends_past_the_model(self):
        lo, hi = self.extent.box(500, 1.0)
        self.assertAlmostEqual(lo[2], 0)
        self.assertAlmostEqual(hi[2], 200)

    def test_crop_without_clouds_is_refused(self):
        bare = Extent({}, np.array([[0., 0, 0], [1, 1, 1]]))
        self.assertEqual(len(bare.box(0, 0)), 2)
        with self.assertRaises(ValueError):
            bare.box(10, 0)


if __name__ == '__main__':
    unittest.main()
