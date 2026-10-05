"""Unit/coverage contracts independent of the nonlinear simulator."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from field_adapter import trajectory_positions,geometry_hash,load_transfer
from stimulation import Pulse,point_source_basis,combine_basis
from simulation import propagated


class Contracts(unittest.TestCase):
    def test_point_source_units(self):
        # 1 A, 1 metre away, sigma=1 S/m -> 1/(4pi) volts.
        phi=point_source_basis([[1e6,0,0]],[[0,0,0]],1)
        self.assertAlmostEqual(phi[0,0],1/(4*np.pi))
        # 2 mA * V/A -> 2 * phi mV at the NEURON boundary.
        self.assertAlmostEqual(phi[0,0]*2,1/(2*np.pi))

    def test_montage_rejects_unbalanced_nan(self):
        for phi,w in (([[1],[2]],[1,0]),([[np.nan],[2]],[-1,1])):
            with self.assertRaises(ValueError):combine_basis(phi,w)

    def test_biphasic_charge(self):
        for gap in (0,.02):
            t,w=Pulse(width_ms=.3,gap_ms=gap).vectors(10)
            self.assertAlmostEqual(float(np.trapezoid(w,t)),0.,places=12)

    def test_polyline_arclength(self):
        xyz=trajectory_positions([[0,0,0],[1,0,0],[1,2,0]],[500,1500,3000])
        np.testing.assert_allclose(xyz,[[.5,0,0],[1,.5,0],[1,2,0]])
        with self.assertRaises(ValueError):trajectory_positions([[0,0,0],[1,0,0]],[1001])

    def test_contact_identity_and_geometry(self):
        s=np.array([0.,1.]);xyz=np.zeros((2,3))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'field.npz'
            np.savez(p,s_um=s,xyz_mm=xyz,geometry_hash=geometry_hash('c',s,xyz),
                     units_coords='mm',units_phi='V/A',contact_ids=[3,7],phi_V_A=[[1,2],[4,8]],metadata_json='{}')
            transfer,_=load_transfer(p,'c',s,xyz,{3:-1,7:1})
            np.testing.assert_allclose(transfer,[3,6])
            with self.assertRaises(ValueError):load_transfer(p,'c',s,xyz,{1:-1,7:1})
            with self.assertRaises(ValueError):load_transfer(p,'c',s,xyz+1,{3:-1,7:1})

    def test_simultaneous_excursions_are_not_propagation(self):
        t=np.arange(10.)
        v=np.full((9,10),-60.);v[:,5:7]=20
        self.assertFalse(propagated(t,v,np.arange(9.),0)[0])
        v[7,:]=-60;v[7,6:8]=20
        self.assertTrue(propagated(t,v,np.arange(9.),0)[0])


if __name__ == '__main__':unittest.main()
