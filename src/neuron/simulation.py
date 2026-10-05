"""Repeatable single-pulse trials and conservative propagated-spike detection."""
import numpy as np
from neuron import h
from stimulation import Pulse


def crossings(t, v, after_ms=0.):
    i = np.flatnonzero((v[:-1] < 0) & (v[1:] >= 0))
    times = t[i] + (-v[i])*(t[i+1]-t[i])/(v[i+1]-v[i])
    return times[times >= after_ms]


def propagated(t, traces, positions, after_ms):
    times = [crossings(t,v,after_ms) for v in traces]
    # Two successive sites on either side of the stimulation region must
    # cross 0 mV in outward order. A local excursion alone is not activation.
    arrivals = []
    for a,b in ((2,1),(6,7)):
        for ta in times[a]:
            later = times[b][times[b] > ta + 1e-6]
            if later.size:
                tb = float(later[0])
                arrivals.append(dict(sites=[a,b],times_ms=[float(ta),tb],
                    velocity_m_s=float(abs(positions[b]-positions[a])/(tb-ta)*.001)))
                break
    return bool(arrivals), arrivals


def trial(axon, transfer_V_A=None, amplitude_mA=0., pulse=Pulse(), dt_ms=.005,
          stop_ms=None, intracellular_nA=None):
    if not np.isfinite(dt_ms) or dt_ms <= 0 or dt_ms > pulse.width_ms/10:
        raise ValueError('dt must be positive and resolve each pulse phase in at least 10 steps')
    if not np.isfinite(amplitude_mA) or amplitude_mA < 0:
        raise ValueError('amplitude must be finite and nonnegative')
    # Allow distant recording sites to be reached when the caller lengthens
    # the axon. These are conservative run-window estimates, not measured CV.
    speed_floor = {'abeta':10., 'adelta':2., 'c':.3}[axon.kind]
    stop_ms = stop_ms or max(40. if axon.kind == 'c' else 15.,
                            pulse.start_ms+2*pulse.width_ms+pulse.gap_ms+
                            axon.length_um*.001/speed_floor+5.)
    tv, waveform = pulse.vectors(stop_ms)
    if transfer_V_A is not None:
        transfer_V_A = np.asarray(transfer_V_A,float)
        if transfer_V_A.shape != axon.s_um.shape or not np.isfinite(transfer_V_A).all():
            raise ValueError('One finite transfer value is required per compartment')
    plays, clamps = [], []
    time_play = h.Vector(tv)
    sites, positions = axon.recording_sites(tuple(np.arange(1,10)/10))
    records = [h.Vector().record(seg._ref_v) for seg in sites]
    clock = h.Vector().record(h._ref_t)
    try:
        for j,seg in enumerate(axon.segments):
            seg.e_extracellular = 0
            if transfer_V_A is not None:
                # V/A * (mA * 1e-3) * 1e3 mV/V: numeric scale is one.
                play = h.Vector(waveform * amplitude_mA * transfer_V_A[j])
                play.play(seg._ref_e_extracellular, time_play, 1)
                plays.append(play)
        if intracellular_nA is not None:
            site = axon.recording_sites((.05,))[0][0]
            stim = h.IClamp(site)
            stim.delay, stim.dur, stim.amp = pulse.start_ms, pulse.width_ms, intracellular_nA
            clamps.append(stim)
        h.CVode().active(0)
        h.secondorder = 0
        h.celsius, h.dt = axon.temperature_C, dt_ms
        h.steps_per_ms = 1/dt_ms
        h.finitialize(axon.rest_mV)
        h.continuerun(stop_ms)
        t, v = np.array(clock), np.array([np.array(r) for r in records])
        if not np.isfinite(v).all():
            raise RuntimeError('Nonfinite membrane voltage')
        active, arrivals = propagated(t,v,positions,pulse.start_ms)
        spontaneous = any(crossings(t,row,0).size != crossings(t,row,pulse.start_ms).size for row in v)
        return dict(t_ms=t,v_mV=v,s_um=positions,activated=active,
                    arrivals=arrivals,spontaneous=spontaneous,
                    peak_mV=float(v.max()),dt_ms=dt_ms)
    finally:
        for play in plays:
            play.play_remove()
        for seg in axon.segments:
            seg.e_extracellular = 0


def threshold(axon, transfer, pulse=Pulse(), dt_ms=.005, max_mA=100., rtol=.01):
    if not np.isfinite(max_mA) or max_mA <= 0 or not 0 < rtol < 1:
        raise ValueError('Invalid threshold limits')
    history = []
    def evaluate(amp):
        result = trial(axon,transfer,amp,pulse,dt_ms)
        if result['spontaneous']:
            raise RuntimeError('Spontaneous spike before stimulation; threshold undefined')
        history.append(dict(amplitude_mA=float(amp),activated=result['activated']))
        return result
    baseline = evaluate(0.)
    if baseline['activated']:
        raise RuntimeError('Zero-current control propagates; threshold undefined')
    low, high = 0., min(.001,max_mA)
    while True:
        result = evaluate(high)
        if result['activated']:
            break
        if high == max_mA:
            return dict(status='no_activation_within_range',threshold_mA=None,
                        tested_max_mA=max_mA,history=history),result
        low, high = high, min(high*2,max_mA)
    while high-low > max(1e-7,rtol*high):
        mid = (low+high)/2
        result = evaluate(mid)
        if result['activated']:
            high = mid
        else:
            low = mid
    # Explicitly retain the activated upper endpoint, never the last midpoint.
    result = evaluate(high)
    return dict(status='bracketed',threshold_mA=high,lower_mA=low,upper_mA=high,
                relative_bracket=(high-low)/high,history=history,
                criterion='outward 0-mV crossings at two separated sites'),result
