"""Regressions for idle montage recomputes and Gmsh worker configuration."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import scs_montage_feature as F


class Feature:
    def __init__(self):
        self.PropertiesList = ['SolutionNpz', 'Status', 'ClampVolts']
        self.SolutionNpz = ''
        self.Status = ''
        self.ClampVolts = 0.

    def addProperty(self, kind, name, *args):
        self.PropertiesList.append(name)
        setattr(self, name, 0 if kind == "App::PropertyInteger" else "")

    def setEditorMode(self, *args):
        pass

    def removeProperty(self, name):
        self.PropertiesList.remove(name)
        delattr(self, name)


class VisualizerThreadsTests(unittest.TestCase):
    def test_idle_recompute_never_loads_reference_solution(self):
        for count in (0, 4, 6, 8):
            fp = Feature()
            F.SCSMontageProxy.ensure_contact_properties(fp, count)
            with patch.object(F.M, 'run_montage') as run:
                F.SCSMontageProxy.execute(None, fp)
                F.SCSMontageViewProxy.updateData(None, fp, 'Status')
                run.assert_not_called()
            self.assertIn('No field selected', fp.Status)

    def test_contact_changes_unbind_previous_solution(self):
        fp = Feature()
        for count in (8, 4, 6, 0):
            fp.SolutionNpz = '/old/solution.npz'
            F.SCSMontageProxy.ensure_contact_properties(fp, count)
            self.assertEqual(F.contact_indices(fp), list(range(1, count + 1)))
            self.assertEqual(fp.SolutionNpz, '')

    def test_bound_montage_passes_selected_contact_count_and_file(self):
        for count in (4, 6, 8, 16):
            fp = Feature()
            F.SCSMontageProxy.ensure_contact_properties(fp, count)
            fp.SolutionNpz = '/selected-lead/solution.npz'
            stats = dict(clamp_v=1., charge_sum_mA=0., tissues={
                t: dict(v_min=0., v_max=0., n_nan=0) for t in F.M.TISSUES})
            proxy = F.SCSMontageProxy.__new__(F.SCSMontageProxy)
            with patch.object(proxy, 'validate_binding'), patch.object(F.M, 'run_montage', return_value=stats) as run:
                proxy.execute(fp)
                self.assertEqual(len(run.call_args.args[0]), count)
                self.assertEqual(run.call_args.kwargs['sol_npz'], fp.SolutionNpz)

    def test_real_gmsh_build_receives_twenty_thread_limit(self):
        import gmsh
        import build_mesh
        seen = []
        generate = gmsh.model.mesh.generate

        def generate_checked(dim):
            seen.append([gmsh.option.getNumber(name) for name in (
                'General.NumThreads', 'Mesh.MaxNumThreads1D',
                'Mesh.MaxNumThreads2D', 'Mesh.MaxNumThreads3D')])
            return generate(dim)

        params = dict(H_MIN=.3, H_MAX=.3, FIELDS=[], GRID=.5, MARGIN=0.)
        bounds = np.array([[0., 0., 0.], [2., 2., 2.]])
        with tempfile.TemporaryDirectory() as out, \
                patch.dict(os.environ, {'OMP_NUM_THREADS': '20'}), \
                patch('mesh_cost.CostModel'), \
                patch('mesh_cost.describe', return_value=('ok', 'test cube')), \
                patch('mesh_geometry.geometry', return_value=({}, bounds)), \
                patch.object(gmsh.model.mesh, 'generate', side_effect=generate_checked):
            try:
                nodes, tets, _, _ = build_mesh.build(params, out_dir=out)
            finally:
                if gmsh.isInitialized():
                    gmsh.finalize()
        # Boundary meshing intentionally stays serial to avoid the Gmsh crash.
        self.assertEqual(seen, [[20., 1., 1., 20.]])
        self.assertGreater(len(tets), 0)
        xyz = nodes[tets]
        volumes = np.abs(np.linalg.det(xyz[:, 1:] - xyz[:, :1])) / 6
        self.assertAlmostEqual(float(volumes.sum()), 8., places=6)


if __name__ == '__main__':
    unittest.main()
