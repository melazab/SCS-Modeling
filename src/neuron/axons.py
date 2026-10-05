"""Representative straight axons; see MODEL_SOURCES.md for scope and provenance.

Lengths/diameters are micrometres. No global hoc geometry, stimulus, or GUI.
Use one model per process: NEURON temperature and mechanism state are global.
"""
import math
import numpy as np
from neuron import h
from build_mechanisms import build

_LOADED = False


class Axon:
    def __init__(self, kind, length_um=10000., refinement=1, temperature_C=37.):
        global _LOADED
        if kind not in ('abeta', 'adelta', 'c'):
            raise ValueError('kind must be abeta, adelta, or c')
        if not np.isfinite(length_um) or length_um < 3000:
            raise ValueError('length must be finite and at least 3000 um')
        if refinement not in (1, 2):
            raise ValueError('refinement must be 1 or 2')
        if not _LOADED:
            h.nrn_load_dll(str(build()))
            h.load_file('stdrun.hoc')
            _LOADED = True
        self.kind = kind
        self.sections, self.nodes = [], []
        self.rest_mV = -80. if kind == 'abeta' else -60.
        if not np.isfinite(temperature_C) or not 20 <= temperature_C <= 40:
            raise ValueError('temperature must be between 20 and 40 C')
        self.temperature_C = float(temperature_C)
        self.diameter_um = {'abeta': 5.7, 'adelta': 3., 'c': .8}[kind]
        h.celsius = self.temperature_C
        if kind == 'abeta':
            self._mrg(length_um, refinement)
        elif kind == 'adelta':
            self._adelta(length_um, refinement)
        else:
            self._c(length_um, refinement)
        self.length_um = sum(s.L for s in self.sections)
        self.segments, positions = [], []
        self.node_positions = []
        offset = 0.
        for sec in self.sections:
            length, diameter = sec.L, sec.diam
            h.pt3dclear(sec=sec)
            h.pt3dadd(offset, 0, 0, diameter, sec=sec)
            h.pt3dadd(offset + length, 0, 0, diameter, sec=sec)
            for seg in sec:
                self.segments.append(seg)
                positions.append(offset + seg.x * sec.L)
            if sec in self.nodes:
                self.node_positions.append(offset + sec.L / 2)
            offset += sec.L
        self.s_um = np.asarray(positions)
        self.node_positions = np.asarray(self.node_positions)
        # Published A-delta/Sundt initialization balances leak against ionic
        # currents at rest. Do it once, then hold e_pas fixed across trials.
        h.finitialize(self.rest_mV)
        h.fcurrent()
        if kind != 'abeta':
            for seg in self.segments:
                if hasattr(seg, 'ttxs'):
                    seg.pas.e = self.rest_mV + seg.ina / seg.pas.g
                elif hasattr(seg, 'nahh'):
                    seg.pas.e = self.rest_mV + (seg.ina + seg.ik) / seg.pas.g

    def _section(self, name, length, diam, nseg=1):
        sec = h.Section(name=f'{self.kind}_{name}')
        sec.L, sec.diam, sec.nseg = length, diam, nseg
        sec.insert('extracellular')
        if self.sections:
            sec.connect(self.sections[-1](1))
        self.sections.append(sec)
        return sec

    def _mrg(self, length, refinement):
        # Exact 5.7 um table entry in upstream MRGaxon.hoc; no extrapolation.
        fiber, axon, node, dx, para2, nl = 5.7, 3.4, 1.9, 500., 35., 80
        count = int(math.ceil((length - 1) / dx)) + 1
        def configure(sec, inner, gap, gpas=None):
            sec.Ra = 70 * (fiber / inner)**2 if gpas is not None else 70
            sec.cm = 2 * inner / fiber if gpas is not None else 2
            if gpas is not None:
                sec.insert('pas')
                sec.g_pas, sec.e_pas = gpas * inner / fiber, -80
            else:
                sec.insert('axnode')
            for seg in sec:
                seg.xraxial[0] = 7000 / (math.pi * ((inner / 2 + gap)**2 - (inner / 2)**2))
                seg.xg[0] = .001 / (2*nl) if gpas is not None else 1e10
                seg.xc[0] = .1 / (2*nl) if gpas is not None else 0
        for i in range(count):
            sec = self._section(f'node{i}', 1., node, 3 if refinement == 1 else 5)
            self.nodes.append(sec)
            configure(sec, node, .002)
            if i == count - 1:
                break
            pieces = [('MYSA', 3., node, .002, .001), ('FLUT', para2, axon, .004, .0001)]
            pieces += [('STIN', (dx-1-6-2*para2)/6, axon, .004, .0001)] * 6
            pieces += [('FLUT', para2, axon, .004, .0001), ('MYSA', 3., node, .002, .001)]
            for j,(name, L, inner, gap, gpas) in enumerate(pieces):
                sec = self._section(f'{name}{i}_{j}', L, fiber, 3 if refinement == 1 else 5)
                configure(sec, inner, gap, gpas)

    def _adelta(self, length, refinement):
        # Central (dorsal-root) axon of Hao et al. 2023 / DRGsims.
        count = int(math.ceil((length-1.5)/151.5)) + 1
        for i in range(count):
            sec = self._section(f'node{i}', 1.5, 3., 5 if refinement == 1 else 11)
            self.nodes.append(sec)
            sec.Ra, sec.cm = 100, .5
            sec.insert('ttxs')
            sec.gbar_ttxs = .6
            # NEURON's default na_ion reversal used by upstream code.
            sec.ena = 50
            sec.insert('pas')
            sec.g_pas, sec.e_pas = 1/5000, -60
            if i < count-1:
                sec = self._section(f'internode{i}', 150, 3., 25 if refinement == 1 else 51)
                sec.Ra, sec.cm = 100, .01
                sec.insert('pas')
                sec.g_pas, sec.e_pas = 1/100000, -60

    def _c(self, length, refinement):
        # Sundt Figure5d peripheral axon, without soma/T-junction mechanisms.
        nseg = int(math.ceil(length / (50/refinement)))
        nseg += 1 - nseg % 2
        sec = self._section('cable', length, .8, nseg)
        sec.Ra, sec.cm = 100, 1
        sec.insert('nahh')
        sec.gnabar_nahh, sec.mshift_nahh, sec.hshift_nahh = .04, -6, 6
        sec.insert('borgkdr')
        sec.gkdrbar_borgkdr = .04
        sec.ena, sec.ek = 50, -90
        sec.insert('pas')
        sec.g_pas, sec.e_pas = 1/10000, -60

    def recording_sites(self, fractions=(.2,.4,.6,.8)):
        if self.nodes:
            idx = [int(np.argmin(abs(self.node_positions-f*self.length_um))) for f in fractions]
            return [self.nodes[i](.5) for i in idx], self.node_positions[idx]
        idx = [int(np.argmin(abs(self.s_um-f*self.length_um))) for f in fractions]
        return [self.segments[i] for i in idx], self.s_um[idx]

    def close(self):
        self.segments.clear()
        self.nodes.clear()
        for sec in self.sections:
            h.delete_section(sec=sec)
        self.sections.clear()
