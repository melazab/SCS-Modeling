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


if __name__ == '__main__':
    unittest.main()
