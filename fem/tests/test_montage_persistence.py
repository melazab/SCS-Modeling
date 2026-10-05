import sys
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import solve_cpu_limit
import scs_montage as M
import scs_montage_feature as F
import scs_waveform as W


class MontagePersistenceTests(unittest.TestCase):
    def test_physical_limit_counts_unique_cores_in_affinity(self):
        def topology(path):
            if path.name == 'physical_package_id':
                return '0'
            cpu = int(path.parent.parent.name[3:])
            return str(cpu // 2 if cpu < 12 else cpu - 6)
        with patch.object(solve_cpu_limit.os, 'sched_getaffinity', return_value=set(range(20))), \
                patch.object(Path, 'read_text', topology):
            self.assertEqual(solve_cpu_limit.physical_core_limit(), 14)
        with patch.object(solve_cpu_limit.os, 'sched_getaffinity', return_value={0, 1, 12}), \
                patch.object(Path, 'read_text', topology):
            self.assertEqual(solve_cpu_limit.physical_core_limit(), 2)

    def test_recompute_notifies_even_when_summary_is_identical(self):
        fp = SimpleNamespace(SolutionNpz='/run/solution.npz', ClampVolts=0.,
                             Contact1_mA=-1., Contact2_mA=1., PlotRevision=0, WaveformMode='Static',
                             PropertiesList=['Contact1_mA','Contact2_mA','PlotRevision'],
                             setEditorMode=lambda *args: None)
        proxy = F.SCSMontageProxy.__new__(F.SCSMontageProxy)
        stats = dict(clamp_v=1., charge_sum_mA=0., tissues={
            t: dict(v_min=-1.,v_max=1.,n_nan=0) for t in M.TISSUES})
        with patch.object(proxy, 'validate_binding'), patch.object(W, 'ensure_properties'), patch.object(M, 'run_montage', return_value=stats):
            proxy.execute(fp)
            status = fp.Status
            fp.Contact1_mA, fp.Contact2_mA = 1., -1.
            proxy.execute(fp)
        self.assertEqual(fp.Status, status)
        self.assertEqual(fp.PlotRevision, 2)

    def test_biphasic_snapshot_keeps_peak_currents_and_scale(self):
        fp = SimpleNamespace(SolutionNpz='/run/solution.npz', ClampVolts=0.,
             Contact1_mA=-2.4, Contact2_mA=2.4, PlotRevision=0,
             WaveformMode='FAST-style biphasic', FrequencyHz=90., PhaseWidthUs=220.,
             InterphaseGapUs=0., WaveformTimeMs=.3,
             PropertiesList=['Contact1_mA','Contact2_mA','PlotRevision'],
             setEditorMode=lambda *args: None)
        proxy=F.SCSMontageProxy.__new__(F.SCSMontageProxy)
        stats=dict(clamp_v=2.,charge_sum_mA=0.,tissues={
            t:dict(v_min=-2.,v_max=2.,n_nan=0) for t in M.TISSUES})
        with patch.object(proxy,'validate_binding'), patch.object(W,'ensure_properties'), \
             patch.object(M,'choose_clamp_v',return_value=2.), \
             patch.object(M,'run_montage',return_value=stats) as plot:
            proxy.execute(fp)
            self.assertEqual(plot.call_args.args[0],[2.4,-2.4])
            self.assertEqual(plot.call_args.kwargs['clamp_v'],2.)
            fp.WaveformTimeMs=5.
            proxy.execute(fp)
            self.assertEqual(plot.call_args.args[0],[0.,0.])
            self.assertEqual(plot.call_args.kwargs['clamp_v'],2.)
        self.assertEqual(fp.Contact1_mA,-2.4)
        self.assertEqual(fp.PlotRevision,2)

    def test_invalid_binding_does_not_plot(self):
        fp = SimpleNamespace(SolutionNpz='/stale/solution.npz')
        proxy = F.SCSMontageProxy.__new__(F.SCSMontageProxy)
        with patch.object(proxy, 'validate_binding', side_effect=RuntimeError('Lead geometry changed')), \
                patch.object(M, 'run_montage') as plot:
            proxy.execute(fp)
            plot.assert_not_called()
        self.assertEqual(fp.Status, 'Lead geometry changed')

    def test_lead_surfaces_come_from_selected_run(self):
        with patch.object(M.os.path, 'isdir', return_value=True), \
                patch.object(M.C, 'lead_stls', return_value=({1:'a.stl',2:'b.stl'},'i.stl')) as lead:
            paths = M.surface_paths('/my/run/solution_bg.npz')
        lead.assert_called_once_with('/my/run/lead')
        self.assertEqual(paths['contact_2'], 'b.stl')
        self.assertEqual(paths['insulator'], 'i.stl')
        self.assertNotIn('contact_3', paths)


if __name__ == '__main__':
    unittest.main()
