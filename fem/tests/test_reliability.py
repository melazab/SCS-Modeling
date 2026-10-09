import ast
from collections import Counter
import contextlib
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.spatial import Delaunay

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import artifacts
import config as C
from field import TetField, structured_export
from mesh_faces import face_neighbours, subset_neighbours, face_chunks
from mesh_preview import tissue_surfaces, outer_surface
from assign_and_solve import dura_leak_report, check_residual, solve_contact
from mesh_cost import describe
from export_results import write_vtu
from mesh_geometry import grid_points


class ReliabilityTests(unittest.TestCase):
    def test_solver_only_change_reuses_mesh_but_invalidates_solution(self):
        original = artifacts.file_hash
        mesh = artifacts.model_signature()
        solve = artifacts.model_signature('solve')
        with patch.object(artifacts, 'file_hash', side_effect=lambda p: 'changed-solver' if str(p).endswith('/solve_lead.py') else original(p)):
            self.assertEqual(mesh, artifacts.model_signature())
            self.assertNotEqual(solve, artifacts.model_signature('solve'))

    def test_cg_continues_without_relaxing_acceptance(self):
        from unittest.mock import Mock
        from scipy.sparse import eye
        ml = Mock()
        ml.solve.side_effect = [np.zeros(2), np.ones(2)]
        x, residual, iterations = solve_contact(ml, eye(2), np.ones(2), 1)
        self.assertEqual(ml.solve.call_count, 2)
        self.assertEqual(residual, 0)
        np.testing.assert_array_equal(ml.solve.call_args.kwargs['x0'], np.zeros(2))

    def test_cg_rejects_after_all_batches_fail(self):
        from unittest.mock import Mock
        from scipy.sparse import eye
        ml = Mock()
        ml.solve.return_value = np.zeros(2)
        with self.assertRaisesRegex(RuntimeError, 'no solution published'):
            solve_contact(ml, eye(2), np.ones(2), 1)
        self.assertEqual(ml.solve.call_count, 4)

    def test_direct_solver_enforces_memory_budget_before_classification(self):
        from assign_and_solve import solve_basis_fields
        with patch('mesh_cost.available_ram_gb', return_value=0.1), patch('assign_and_solve.classify_tets') as classify:
            with self.assertRaisesRegex(RuntimeError, 'coarser mesh'):
                solve_basis_fields(np.zeros((4, 3)), np.array([[0, 1, 2, 3]]))
            classify.assert_not_called()

    def test_direct_mesher_enforces_memory_budget_before_geometry(self):
        from build_mesh import build
        with tempfile.TemporaryDirectory() as tmp, patch('artifacts.valid', return_value=False), patch('mesh_cost.CostModel') as model, patch('mesh_cost.describe', return_value=('danger', 'too large')), patch('mesh_geometry.geometry') as geometry:
            model.return_value.estimate.return_value = {}
            with self.assertRaisesRegex(RuntimeError, 'refused before build'):
                build(out_dir=tmp)
            geometry.assert_not_called()

    def test_unconverged_and_nonfinite_fields_are_rejected(self):
        from scipy.sparse import csr_matrix
        a = csr_matrix([[2., -1.], [-1., 2.]])
        rhs = np.array([1., 0.])
        self.assertLess(check_residual(a, np.linalg.solve(a.toarray(), rhs), rhs, 1, 2), 1e-12)
        for bad in (np.zeros(2), np.array([np.nan, 1.])):
            with self.assertRaisesRegex(RuntimeError, 'no solution published'):
                check_residual(a, bad, rhs, 1, 400)

    def test_face_neighbours_match_independent_face_count(self):
        nodes = np.random.default_rng(31).random((80, 3))
        tets = Delaunay(nodes).simplices
        expected = Counter()
        for tet in tets:
            for omit in range(4):
                expected[tuple(sorted(np.delete(tet, omit)))] += 1
        nbr = face_neighbours(tets, threads=3)
        found = Counter()
        lone = set()
        for f, owner, other in face_chunks(tets, nbr, chunk=13):
            for face, who, across in zip(f, owner, other):
                key = tuple(sorted(face))
                found[key] += 1
                self.assertTrue(set(face).issubset(tets[who]))
                if across < 0:
                    lone.add(key)
                else:
                    self.assertTrue(set(face).issubset(tets[across]))
                    self.assertIn(who, nbr[across])
        self.assertEqual(set(found.values()), {1})
        self.assertEqual(set(found), set(expected))
        self.assertEqual(lone, {f for f, n in expected.items() if n == 1})
        keep = np.arange(len(tets)) % 3 != 0
        np.testing.assert_array_equal(subset_neighbours(nbr, keep), face_neighbours(tets[keep]))

    def test_face_neighbours_reject_non_manifold_mesh(self):
        with self.assertRaisesRegex(ValueError, 'Non-manifold'):
            face_neighbours(np.array([[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 2, 5]]))

    def test_memory_guard_does_not_count_forked_shared_pages(self):
        # Parallel classification forks workers that read the parent's arrays.
        # Summed RSS counted those shared pages once per worker and killed solves.
        script = """
import sys, multiprocessing, numpy as np
sys.path.insert(0, %r)
from solve_resources import Resources
BIG = None
def work(i):
    import time
    time.sleep(1)  # outlive a 0.25 s guard sample, like real classification chunks
    return float(BIG[i::97].sum())
with Resources(1, 1.0) as r:
    BIG = np.ones(int(0.3e9 / 8))
    with multiprocessing.get_context('fork').Pool(8) as pool:
        pool.map(work, range(8))
print('peak %%.2f' %% r.peak_gb)
""" % str(Path(__file__).resolve().parents[1] / 'scripts')
        import subprocess
        out = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True, timeout=120)
        self.assertEqual(out.returncode, 0, out.stdout + out.stderr)
        self.assertLess(float(out.stdout.split('peak')[-1]), 1.0)

    def test_tissue_interfaces_and_leak_area(self):
        nodes = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1]])
        tets = np.array([[0, 1, 2, 3], [0, 2, 1, 4]])
        order = ['epidural', 'csf', 'white', 'grey', 'dura', 'insulator']
        lab = np.array([0, 1])
        f, ids = tissue_surfaces(nodes, tets, lab)
        self.assertEqual(len(f), 8)  # six external, two sides of interface
        self.assertEqual(Counter(ids), {0: 4, 1: 4})
        f, ids = outer_surface(nodes, tets, lab)
        self.assertEqual(Counter(ids), {0: 3, 1: 3})
        self.assertAlmostEqual(dura_leak_report(nodes, tets, lab, order)['leak_mm2'], 0.5)

    def test_small_mesh_sampling_and_outside(self):
        n = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        f = TetField(n, [[0, 1, 2, 3]], n.sum(1))
        v = f([[.1, .2, .3], [3., 3, 3]])
        self.assertAlmostEqual(v[0, 0], .6)
        self.assertTrue(np.isnan(v[0, 1]))
        self.assertEqual(f(np.empty((0, 3))).shape, (1, 0))

    def test_four_contact_exports(self):
        import meshio
        n = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        phi = np.arange(1, 5)[:, None] * n.sum(1)
        d = dict(nodes=n, tets=np.array([[0, 1, 2, 3]]), phi=phi, V=phi[0],
                 label=np.array([0]), sigma=np.array([1.]), order=np.array(['white']))
        with tempfile.TemporaryDirectory() as tmp:
            sol, grid, vtu = [os.path.join(tmp, x) for x in ('sol.npz', 'grid.npz', 'sol.vtu')]
            np.savez(sol, **d)
            structured_export(sol, grid, [.1, .1, .1], .1, (2, 2, 2))
            with np.load(grid) as g:
                self.assertEqual(g['phi'].shape, (4, 2, 2, 2))
                np.testing.assert_allclose(g['phi'][:, 0, 0, 0], [.3, .6, .9, 1.2])
            write_vtu(d, vtu)
            self.assertEqual(len(meshio.read(vtu).point_data), 5)

    def test_artifact_corruption_and_input_changes_are_cache_misses(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, 'solution.npz')
            artifacts.atomic_npz(p, V=[1., 2.])
            self.assertFalse(artifacts.valid(p, 'model-a'))
            artifacts.publish(p, 'model-a')
            self.assertTrue(artifacts.valid(p, 'model-a'))
            self.assertFalse(artifacts.valid(p, 'model-b'))
            artifacts.atomic_npz(p, V=[4., 5.])
            self.assertFalse(artifacts.valid(p, 'model-a'))

    def test_signatures_track_material_and_geometry_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            tissue = Path(tmp) / 'tissue.yaml'
            anatomy = Path(tmp) / 'anatomy.stl'
            tissue.write_text('sigma: 1')
            anatomy.write_text('geometry-a')
            with patch.object(C, 'TISSUE_MAP_YAML', str(tissue)), patch.object(C, 'TISSUE_BODIES', {'x': [str(anatomy)]}):
                old = artifacts.model_signature()
                tissue.write_text('sigma: 2')
                changed = artifacts.model_signature()
                self.assertNotEqual(old, changed)
                anatomy.write_text('geometry-b')
                self.assertNotEqual(changed, artifacts.model_signature())

    def test_background_outputs_are_separate(self):
        with patch.dict(os.environ, {'SCS_BACKGROUND': '1'}):
            for name in ['solution.npz', 'field_grid.npz', 'mesh_tagged.msh', 'voltage_slices.png']:
                self.assertIn('_bg.', artifacts.result_path(name))
        with patch.dict(os.environ, {'SCS_BACKGROUND': '0'}):
            self.assertEqual(artifacts.result_path('solution.npz'), os.path.join(C.OUT, 'solution.npz'))

    def test_analytic_gate_does_not_mask_a_failed_resolution(self):
        source = ast.parse((Path(__file__).parents[1] / 'scripts/verify_solver.py').read_text())
        exit_expr = source.body[-1].body[-1].value.args[0]
        for codes in ([0, 1, 0], [1, 0, 1], [0, 0, 1]):
            code = compile(ast.Expression(exit_expr), '<verification exit>', 'eval')
            self.assertEqual(eval(code, {'codes': codes}), 1)

    def test_memory_gate_uses_upper_estimate(self):
        estimate = dict(n_tets=47_500_000, ram_gb_lo=19.5, ram_gb_hi=34.25, mesh_seconds=2000)
        self.assertEqual(describe(estimate, ram_gb=24)[0], 'danger')

    def test_grid_chunks_preserve_structured_order(self):
        lo, n, step = np.array([1., 2., 3.]), (3, 4, 5), .25
        full = np.stack(np.meshgrid(*[lo[i] + step * np.arange(n[i]) for i in range(3)], indexing='ij'), -1).reshape(-1, 3)
        chunked = np.concatenate([grid_points(lo, n, step, s, min(s + 7, 60)) for s in range(0, 60, 7)])
        np.testing.assert_array_equal(full, chunked)


if __name__ == '__main__':
    unittest.main()
