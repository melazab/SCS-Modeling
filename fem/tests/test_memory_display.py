"""Memory-gate wording and release of closed-document FEM caches."""
from pathlib import Path
import sys
import unittest
from unittest.mock import Mock, patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import scs_montage as M
import mesh_cost

class MemoryDisplayTests(unittest.TestCase):
    def test_closed_solution_released_open_solution_kept(self):
        closed,opened=Mock(),Mock()
        cache={('/tmp/closed.npz',1):(object(),closed,object()),('/tmp/open.npz',2):(object(),opened,object())}
        with patch.object(M,'_CACHE',cache):
            self.assertEqual(M.release_unused_basis(['/tmp/open.npz']),1)
            self.assertEqual(list(cache),[('/tmp/open.npz',2)])
            self.assertEqual(M.release_unused_basis(['/tmp/open.npz']),0)
        closed.close.assert_called_once()
        opened.close.assert_not_called()

    def test_available_ram_is_not_the_usable_budget(self):
        est=dict(n_tets=17.5e6,mesh_seconds=780,ram_gb_lo=8,ram_gb_hi=19)
        severity,text=mesh_cost.describe(est,ram_gb=20)
        self.assertEqual(severity,'danger')
        self.assertIn('16.0 GB mesh budget',text)
        self.assertIn('20.0 GB available',text)
        self.assertEqual(mesh_cost.describe(est,ram_gb=25)[0],'warn')

if __name__=='__main__':unittest.main()
