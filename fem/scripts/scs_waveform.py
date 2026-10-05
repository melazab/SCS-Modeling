"""Idealized FAST-style current waveform; no electrode/tissue capacitance.

90 Hz and symmetric active recharge follow Metzger 2021. The 220 us default
is the dorsal setting in Van Acker & Kim 2023, used identically for the matched
comparison. Zero interphase gap is an explicit assumption, not a device spec.
"""
import math


def validate(frequency_hz, width_us, gap_us):
    if not all(math.isfinite(v) for v in (frequency_hz, width_us, gap_us)):
        raise ValueError('Waveform parameters must be finite')
    if frequency_hz <= 0 or width_us <= 0 or gap_us < 0:
        raise ValueError('Frequency and phase width must be positive; gap must be nonnegative')
    period_ms = 1000. / frequency_hz
    if (2 * width_us + gap_us) / 1000. >= period_ms:
        raise ValueError('Both phases and the gap must fit inside one period')
    return period_ms


def phase(time_ms, frequency_hz=90., width_us=220., gap_us=0.):
    period = validate(frequency_hz, width_us, gap_us)
    if not math.isfinite(time_ms) or time_ms < 0:
        raise ValueError('Time must be finite and nonnegative')
    t = time_ms % period
    width, gap = width_us / 1000., gap_us / 1000.
    if t < width:
        return 1., 'First phase'
    if t < width + gap:
        return 0., 'Interphase gap'
    if t < 2 * width + gap:
        return -1., 'Reverse phase'
    return 0., 'Interpulse interval'


def totals(currents):
    if not all(math.isfinite(v) for v in currents):
        raise ValueError('Currents must be finite')
    return sum(max(v, 0.) for v in currents), -sum(min(v, 0.) for v in currents)


def require_balanced(currents):
    source, sink = totals(currents)
    if source <= 0 or sink <= 0:
        raise ValueError('Select at least one source and one sink')
    if abs(source-sink) > 1e-8 * max(1., source, sink):
        raise ValueError('Source and sink totals must match; use Scale to total')
    return source


def normalize(currents, total_ma):
    source, sink = totals(currents)
    if not math.isfinite(total_ma) or total_ma <= 0 or source <= 0 or sink <= 0:
        raise ValueError('Choose a positive total and at least one source and sink')
    return [v * total_ma / (source if v > 0 else sink) for v in currents]


def ensure_properties(obj):
    for name, value, description in (
        ('FrequencyHz', 90., 'Pulse repetition frequency (Hz)'),
        ('PhaseWidthUs', 220., 'Duration of EACH equal rectangular phase (microseconds)'),
        ('InterphaseGapUs', 0., 'Assumed interphase gap; zero is not a confirmed device specification'),
        ('WaveformTimeMs', 0., 'Time within the repeating waveform (milliseconds)'),
    ):
        if name not in obj.PropertiesList:
            obj.addProperty('App::PropertyFloat', name, 'Waveform', description)
            setattr(obj, name, value)
    if 'WaveformMode' not in obj.PropertiesList:
        obj.addProperty('App::PropertyEnumeration', 'WaveformMode', 'Waveform')
        obj.WaveformMode = ['Static', 'FAST-style biphasic']
        obj.WaveformMode = 'Static'
