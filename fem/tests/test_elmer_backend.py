import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from element_batches import assemble, action
from solve_resources import Resources
import elmer_backend as E
import fibres


def fixture():
    nodes = np.array(np.meshgrid(*([np.linspace(0, 1, 5)] * 3), indexing='ij')).reshape(3, -1).T
    nodes += np.random.default_rng(15).normal(0, .001, nodes.shape)
    tets = Delaunay(nodes).simplices
    labels = np.clip((nodes[tets].mean(1)[:, 0] * 4).astype(int), 0, 3)
    return nodes, tets, labels


class ElementTests(unittest.TestCase):
    def test_assembly_batches_and_matrix_free_action(self):
        nodes, tets, _ = fixture()
        sigma = np.linspace(.05, 2., len(tets))
        a, vol, _ = assemble(nodes, tets, sigma, batch=17)
        b, v2, _ = assemble(nodes, tets, sigma, batch=113)
        np.testing.assert_allclose(a.toarray(), b.toarray(), rtol=1e-13, atol=1e-15)
        np.testing.assert_allclose(vol, v2)
        values = np.random.default_rng(9).normal(size=len(nodes))
        np.testing.assert_allclose(action(nodes, tets, sigma, values, batch=13), a @ values, rtol=1e-12, atol=1e-15)
        roundoff = 20 * np.finfo(float).eps * float(abs(a).sum(axis=1).max())
        np.testing.assert_allclose(a @ np.ones(len(nodes)), 0, atol=roundoff)

    def test_known_tetrahedron_and_reversed_orientation(self):
        nodes = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        gradients = np.array([[-1., -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]]) * 1000
        expected = 2e-9 / 6 * (gradients @ gradients.T)
        for tet in ([[0, 1, 2, 3]], [[0, 2, 1, 3]]):
            a, vol, _ = assemble(nodes, np.array(tet), np.array([2.]), batch=1)
            np.testing.assert_allclose(a.toarray(), expected)
            np.testing.assert_allclose(vol, [1e-9 / 6])

    def test_anisotropic_element_matches_explicit_tensor(self):
        """sigma I + a a^T, assembled through the fibre term, equals vol G S G^T."""
        nodes, tets, _ = fixture()
        rng = np.random.default_rng(3)
        sigma = rng.uniform(.05, 1., len(tets))
        table = rng.normal(size=(5, 3))
        fibre_bin = rng.integers(-1, 5, len(tets)).astype(np.int16)
        a, vol, _ = assemble(nodes, tets, sigma, batch=19, fibre=(fibre_bin, table))
        ref = np.zeros((len(nodes), len(nodes)))
        for e, tet in enumerate(tets):
            p = nodes[tet] * 1e-3
            m = np.column_stack([np.ones(4), p])
            g = np.linalg.inv(m)[1:].T                       # gradients of the 4 shape functions
            s = sigma[e] * np.eye(3)
            if fibre_bin[e] >= 0:
                s = s + np.outer(table[fibre_bin[e]], table[fibre_bin[e]])
            ref[np.ix_(tet, tet)] += abs(np.linalg.det(m)) / 6 * g @ s @ g.T
        np.testing.assert_allclose(a.toarray(), ref, rtol=1e-10, atol=1e-12 * abs(ref).max())
        values = rng.normal(size=len(nodes))
        np.testing.assert_allclose(action(nodes, tets, sigma, values, batch=7, fibre=(fibre_bin, table)),
                                   a @ values, rtol=1e-12, atol=1e-15)

    def test_anisotropic_linear_field_is_exact(self):
        """Analytic: uniform current J0 through a homogeneous anisotropic block.
        V = -(S^-1 J0).x, and A V must equal the boundary flux load J0.n."""
        nodes, tets, _ = fixture()
        d = np.array([.3, -.5, .8]); d /= np.linalg.norm(d)
        st, sl = .083, .6
        S = st * np.eye(3) + (sl - st) * np.outer(d, d)
        J0 = np.array([.4, -.6, 1.])                 # mostly along the fibres
        V = -(nodes * 1e-3) @ np.linalg.solve(S, J0)
        fibre = (np.zeros(len(tets), np.int16), np.sqrt(sl - st) * d[None])
        a, _, _ = assemble(nodes, tets, np.full(len(tets), st), fibre=fibre)
        from assign_and_solve import boundary_nodes
        faces = boundary_nodes(tets)
        p = nodes[faces] * 1e-3
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        n *= np.sign(np.einsum('ij,ij->i', n, p.mean(1) - 0.5e-3))[:, None]   # outward
        load = np.zeros(len(nodes))
        np.add.at(load, faces.ravel(), np.repeat(-(n @ J0) / 2 / 3, 3))      # inflow = -J0.n
        np.testing.assert_allclose(a @ V, load, atol=1e-9 * abs(load).max())
        # and the isotropic transverse value alone would be wrong
        b, _, _ = assemble(nodes, tets, np.full(len(tets), st))
        self.assertGreater(abs(b @ V - load).max(), .1 * abs(load).max())

    def test_direction_bins_stay_within_tolerance(self):
        t = np.radians(np.linspace(-20, 23, 500))
        tangents = np.column_stack([np.zeros_like(t), -np.sin(t), np.cos(t)])
        sample_bin, dirs = fibres.direction_bins(tangents)
        worst = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', tangents, dirs[sample_bin]), -1, 1))).max()
        self.assertLessEqual(worst, fibres.BIN_DEG)
        self.assertTrue(np.all(np.diff(sample_bin) >= 0))
        self.assertGreater(len(dirs), 43 / (2 * fibres.BIN_DEG))

    def test_pieces_not_touching_the_lead_are_dropped_within_limit(self):
        cube = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
        speck = cube * .1 + 5                                 # 1/1000 of the volume
        nodes = np.vstack([cube, speck, [[9., 9, 9]]])        # plus an unused node
        tets = np.array([[0, 1, 2, 3], [4, 5, 6, 7]])
        keep, report = E.drop_islands(nodes, tets, np.array([0, 1]), ['contact1', 'blood'])
        np.testing.assert_array_equal(keep, [True, False])
        self.assertEqual((report['components'], report['kept_components'], report['dropped_tets']), (2, 1, 1))
        self.assertAlmostEqual(report['dropped_by_tissue']['blood'], 1 / 6000, places=6)
        big = np.vstack([cube, cube + 5])                     # an equal-sized piece
        with self.assertRaisesRegex(RuntimeError, 'not connected to the lead'):
            E.drop_islands(big, tets, np.array([0, 1]), ['contact1', 'white'])

    def test_disconnected_domain_rejected(self):
        with self.assertRaisesRegex(RuntimeError, 'disconnected'):
            E.connectivity(np.array([[0, 1, 2, 3], [4, 5, 6, 7]]), 8)

    def test_resource_budget_and_cpu_validation(self):
        with patch('solve_resources.available_gb', return_value=10.):
            self.assertEqual(Resources(1, 0).memory_gb, 8.)
            self.assertEqual(Resources(1, 20).memory_gb, 8.)
            self.assertEqual(Resources(1, 3).memory_gb, 3.)
            with self.assertRaises(ValueError):
                Resources(1, -1)
            with self.assertRaises(ValueError):
                Resources(100000, 1)
            for cpus, memory in ((0, 1), (1, float('nan')), (1, float('inf'))):
                with self.assertRaises(ValueError):
                    Resources(cpus, memory)

    def test_memory_monitor_stops_over_budget_worker(self):
        code = "from solve_resources import Resources; import time\nwith Resources(1,.001): time.sleep(5)"
        env = dict(os.environ, PYTHONPATH=str(Path(E.__file__).parent))
        result = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 2)
        self.assertIn('Memory budget exceeded', result.stdout)


@unittest.skipUnless((E.ELMER_BIN / 'ElmerSolver_mpi').exists(), 'Elmer installation required')
class ElmerIntegrationTests(unittest.TestCase):
    def test_four_contact_mpi_solution_matches_direct_reference(self):
        nodes, tets, labels = fixture()
        order = ['contact1', 'contact2', 'contact3', 'contact4']
        sigma = np.ones(len(tets))
        a, vol, _ = assemble(nodes, tets, sigma, batch=27)
        with tempfile.TemporaryDirectory() as directory, Resources(2, 2) as resources:
            result = E.solve_prepared(nodes, tets, labels, sigma, order, [1, 2, 3, 4], directory, resources)
            # Independent RHS from exterior face areas and each contact volume.
            from assign_and_solve import boundary_nodes
            faces = boundary_nodes(tets)
            p = nodes[faces] * 1e-3
            area = np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1) / 2
            ret = np.zeros(len(nodes)); np.add.at(ret, faces.ravel(), np.repeat(area / 3, 3)); ret /= ret.sum()
            pin = result['pinned_node']
            free = np.arange(len(nodes)) != pin
            for i in range(4):
                src = np.zeros(len(nodes)); mask = labels == i
                np.add.at(src, tets[mask].ravel(), np.repeat(vol[mask] / (4 * vol[mask].sum()), 4))
                ref = np.zeros(len(nodes))
                ref[free] = spsolve(a[free][:, free], (src - ret)[free])
                np.testing.assert_allclose(result['phi'][i], ref, rtol=1e-5, atol=1e-5)
            self.assertLess(max(result['relative_residuals']), 1e-6)
            # Gauge is a partitionable point boundary, not rank-local node IDs.
            sif = (Path(directory) / 'contact1/case.sif').read_text()
            self.assertNotIn('Target Nodes', sif)

    def test_anisotropic_bodies_match_direct_reference(self):
        """Each fibre bin is its own Elmer body with its own 3x3 tensor; the
        result must match a direct solve of the Python assembly."""
        nodes, tets, labels = fixture()
        order = ['contact1', 'contact2', 'white', 'grey']
        labels = np.where(labels == 3, 1, labels)                # contacts at both ends
        labels = np.where((labels == 1) & (nodes[tets].mean(1)[:, 0] < .6), 2, labels)
        sigma = np.where(labels == 2, .083, 1.)
        rng = np.random.default_rng(5)
        dirs = rng.normal(size=(3, 3)); dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        table = np.sqrt(.6 - .083) * dirs
        fibre_bin = np.where(labels == 2, rng.integers(0, 3, len(tets)), -1).astype(np.int16)
        a, vol, _ = assemble(nodes, tets, sigma, fibre=(fibre_bin, table))
        with tempfile.TemporaryDirectory() as directory, Resources(2, 2) as resources:
            result = E.solve_prepared(nodes, tets, labels, sigma, order, [1, 2], directory, resources,
                                      fibre=(fibre_bin, table))
            sif = (Path(directory) / 'contact1/case.sif').read_text()
            self.assertEqual(sif.count('Electric Conductivity(3,3)'), len(np.unique(fibre_bin[fibre_bin >= 0])))
            from assign_and_solve import boundary_nodes
            faces = boundary_nodes(tets)
            p = nodes[faces] * 1e-3
            area = np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1) / 2
            ret = np.zeros(len(nodes)); np.add.at(ret, faces.ravel(), np.repeat(area / 3, 3)); ret /= ret.sum()
            pin = result['pinned_node']
            free = np.arange(len(nodes)) != pin
            for i in range(2):
                src = np.zeros(len(nodes)); mask = labels == i
                np.add.at(src, tets[mask].ravel(), np.repeat(vol[mask] / (4 * vol[mask].sum()), 4))
                ref = np.zeros(len(nodes))
                ref[free] = spsolve(a[free][:, free], (src - ret)[free])
                np.testing.assert_allclose(result['phi'][i], ref, rtol=1e-5, atol=1e-5 * abs(ref).max())
            np.testing.assert_array_equal(result['fibre_bin'], fibre_bin)


if __name__ == '__main__':
    unittest.main()
