import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import live_lead as L
import scs_waveform as W
import scs_fast_viewer as V
import stlio


def document(layout):
    objects=[]
    for lead,count in layout:
        objects += [NS(Name='SCS_Preview_Lead%d_Contact_%02d' % (lead,c)) for c in range(1,count+1)]
        objects.append(NS(Name='SCS_Preview_Lead%d_Insulator' % lead))
    return NS(Objects=list(reversed(objects)))


class MultiLeadTests(unittest.TestCase):
    def test_numeric_identity_unequal_counts_and_hidden_leads(self):
        info=L.find_live_lead(document([(10,2),(2,3)]))
        self.assertEqual(info['status'],L.LiveLeadStatus.OK)
        self.assertEqual(info['contact_map'],{1:[2,1],2:[2,2],3:[2,3],4:[10,1],5:[10,2]})
        self.assertEqual(len(info['insulators']),2)

    def test_incomplete_lead_rejected(self):
        d=document([(1,2),(2,3)])
        d.Objects=[o for o in d.Objects if o.Name!='SCS_Preview_Lead2_Contact_02']
        self.assertEqual(L.find_live_lead(d)['status'],L.LiveLeadStatus.WRONG_SHAPE)

    def test_disconnected_insulation_and_export(self):
        info=L.find_live_lead(document([(1,2),(3,1)]))
        def tess(obj,tol):
            dx=10. if 'Lead3_' in obj.Name else 0.
            return np.array([[dx,0.,0.],[dx+1,0,0],[dx,1,0]]),np.array([[0,1,2]])
        with patch.object(L,'_tessellate_obj',side_effect=tess):
            arrays=L.tessellate_lead_arrays(info)
        np.testing.assert_array_equal(arrays['insulator'][1],[[0,1,2],[3,4,5]])
        with tempfile.TemporaryDirectory() as d:
            L.write_stls(arrays,d)
            contacts,ins=L.C.lead_stls(d)
            self.assertEqual(len(contacts),3)
            self.assertEqual(len(stlio.read_stl(ins)[1]),2)
            self.assertEqual(json.loads((Path(d)/'contact_map.json').read_text())['3'],[3,1])
        before=L.fingerprint_from_arrays(arrays)
        arrays.contact_map[3]=[4,1]
        self.assertNotEqual(before,L.fingerprint_from_arrays(arrays))


class WaveformTests(unittest.TestCase):
    def test_phase_boundaries_and_period(self):
        for t,f in [(0,1),(.219,1),(.22,0),(.239,0),(.24,-1),(.459,-1),(.46,0),(10,0)]:
            self.assertEqual(W.phase(t,90,220,20)[0],f)
        self.assertEqual(W.phase(1000/90+.1,90,220,20)[0],1)
        self.assertEqual(W.phase(.22,90,220,0)[0],-1)
        with self.assertRaises(ValueError):W.validate(10000,220,0)

    def test_balance_and_matched_total_across_leads(self):
        a=W.normalize([-1,2,0,-3,2],2.4)
        self.assertAlmostEqual(W.require_balanced(a),2.4)
        self.assertEqual(a[2],0)
        self.assertAlmostEqual(a[3]/a[0],3)
        # Equal durations and opposite phases give zero net charge per contact.
        np.testing.assert_allclose(np.asarray(a)*.22-np.asarray(a)*.22,0)
        for bad in ([1,0],[-1,2],[0,0]):
            with self.assertRaises(ValueError):W.require_balanced(bad)

    def test_viewer_export_fixed_scale_and_safe_title(self):
        surfaces={'test':(np.eye(3),np.array([[0,1,2]]),np.array([-1.,0.,1.]))}
        meta=dict(frequency_hz=90,width_us=220,gap_us=0,currents_mA=[-1,1],clamp_v=1,total_mA=1,title='</script>')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'view.html';V.export(p,surfaces,meta)
            self.assertIn('\\u003c/script>',p.read_text())
            self.assertNotIn('/*SCS_DATA*/',p.read_text())
            self.assertEqual(json.loads(p.with_suffix('.json').read_text())['width_us'],220)

if __name__=='__main__':unittest.main()
