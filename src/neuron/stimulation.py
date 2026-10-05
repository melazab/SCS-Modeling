"""Extracellular transfer fields in V/A; waveform current in mA, time in ms."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Pulse:
    start_ms: float = 5.
    width_ms: float = .3
    gap_ms: float = .02

    def vectors(self, stop_ms):
        if not all(np.isfinite([self.start_ms, self.width_ms, self.gap_ms, stop_ms])):
            raise ValueError('Pulse times must be finite')
        if self.start_ms <= 0 or self.width_ms <= 0 or self.gap_ms < 0:
            raise ValueError('Invalid pulse timing')
        a = self.start_ms
        b = a + self.width_ms
        c = b + self.gap_ms
        d = c + self.width_ms
        if stop_ms <= d:
            raise ValueError('Simulation must extend beyond the pulse')
        # Duplicate breakpoints preserve square steps under Vector.play(...,1).
        return np.array([0,a,a,b,b,c,c,d,d,stop_ms]), np.array([0,0,1,1,0,0,-1,-1,0,0.])


def point_source_basis(xyz_um, electrodes_um, sigma_S_m=.2):
    xyz, electrodes = np.asarray(xyz_um,float), np.asarray(electrodes_um,float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or electrodes.ndim != 2 or electrodes.shape[1] != 3:
        raise ValueError('Coordinates must have shape (N,3)')
    if not np.isfinite(xyz).all() or not np.isfinite(electrodes).all() or not np.isfinite(sigma_S_m) or sigma_S_m <= 0:
        raise ValueError('Nonfinite coordinates or invalid conductivity')
    r_m = np.linalg.norm(xyz[None,:,:] - electrodes[:,None,:],axis=2)*1e-6
    if np.any(r_m == 0):
        raise ValueError('Point source coincides with a compartment')
    return 1/(4*np.pi*sigma_S_m*r_m)


def combine_basis(phi_V_A, weights):
    phi, weights = np.asarray(phi_V_A,float), np.asarray(weights,float)
    if phi.ndim != 2 or weights.shape != (phi.shape[0],):
        raise ValueError('Expected one weight per contact basis field')
    if not np.isfinite(phi).all() or not np.isfinite(weights).all() or not np.any(weights):
        raise ValueError('Fields must be finite and montage nonzero')
    if abs(weights.sum()) > 1e-10*max(1.,abs(weights).sum()):
        raise ValueError('This workflow requires a spatially balanced montage')
    return weights @ phi
