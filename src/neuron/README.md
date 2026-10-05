# NEURON — three representative sensory axons

The implemented software milestone uses **Aβ, Aδ and C fibres** with a common
extracellular-field and threshold interface. It uses cached **Elmer** fields;
it does not require Ansys or rerun FEM when pulse amplitude changes.

## Current poster scope — September 30, 2026

The immediate deliverable is the matched dorsal/ventral **FEM field comparison**,
not neuronal recruitment. The axon implementation remains available for later
work; its software validation does not establish activation or analgesia in the
poster models. Spatial material verification and mesh convergence are deferred
until after the poster. See [the outline](../../docs/poster_outline.md) and
[FEM roadmap](../../fem/TODO.md).

## Implemented

- Pinned channel mechanisms and provenance for MRG Aβ, Hao/Jaffe Aδ and Sundt C
  representative axons. See [MODEL_SOURCES.md](MODEL_SOURCES.md) for adaptations.
- Straight multicompartment axons, intracellular propagation and repeated trials.
- Analytic bipolar extracellular fields, symmetric charge-balanced biphasic
  pulses, zero-current/uniform-field controls, propagated-spike detection and
  bracketed threshold search. A failure to find activation returns a null
  threshold with the tested amplitudes, not a fabricated threshold.
- Strict FEM sampling at every compartment, with mm/µm and V/A/mA/mV conversions,
  actual contact IDs, provenance checks and optional expected-tissue checks.
- Time/spatial threshold-refinement runner, numerical reports and PNG/PDF traces.

This is an **engineering milestone**, not physiological validation or population
recruitment data. Anatomical Aδ/C trajectories, published response comparisons,
recovery cycles and waveform-train studies remain next steps. The old HH soma
smoke test is an environment check only; it is not the C-fibre model.

## Setup

```bash
python3 -m venv src/neuron/.venv
src/neuron/.venv/bin/pip install -r src/neuron/requirements.txt
src/neuron/.venv/bin/python src/neuron/smoke_test.py
src/neuron/.venv/bin/python src/neuron/build_mechanisms.py
```

A C/C++ compiler is required. Builds/logs live under `src/neuron/build/`, not in
FreeCAD's environment. Channel source files are unmodified and hash-checked.
Use one axon per process because NEURON temperature is global; parallel sweeps
should use independent processes rather than Python threads sharing NEURON.

## Reproduce the analytic validation

```bash
src/neuron/.venv/bin/python -m unittest discover -s src/neuron/tests -v
src/neuron/.venv/bin/python src/neuron/run_validation.py
```

The runner checks each fibre at the base discretization, half the time step,
and finer spatial discretization. The nine independent cases run in separate
NEURON processes; the default is half the detected logical CPUs, capped at nine.
Use `--jobs 1` for sequential execution or `--jobs N` to choose a limit. The
runner requires less than 5% threshold change (an engineering acceptance
criterion). Outputs are under `out/validation/`:
`summary.json`, per-case `report.json` and traces, `validation.png`, and
`validation.pdf`. A failed case exits nonzero and does not publish a new summary.

One case, or a different pulse width:

```bash
src/neuron/.venv/bin/python src/neuron/run_axon.py \
  --kind c --width-ms 0.3 --temperature 37 --out src/neuron/out/c_example
```

All amplitudes are mA. The analytic case uses two point electrodes in a
homogeneous 0.2 S/m medium, 0.5 mm from a straight axon. These thresholds are
**not SCS thresholds**. The current applies in opposite directions on the two
contacts, then reverses for an equal-duration second phase with a 0.02 ms gap.

## Connect a completed Elmer solve

Provide a JSON polyline with `units: "mm"` and `points: [[x,y,z], ...]` in the
FreeCAD document frame. `allowed_tissues` optionally restricts sampled tissue
labels; `anatomically_validated` records review status. The polyline must cover
the complete constructed axon, whose length rounds up to whole internodes.

```bash
src/neuron/.venv/bin/python src/neuron/run_axon.py \
  --kind abeta --out src/neuron/out/fem_fixture/abeta \
  --solution fem/out/lead_runs/b047c359f91480fea3e651e2/solution.npz \
  --trajectory src/neuron/examples/fem_fixture.json \
  --weights '{"1":-1,"3":1}'
```

Repeat with `--kind adelta` or `--kind c` and a different output directory.
Contact IDs are explicit: current on contact *i* is amplitude times its weight.
Weights must sum to zero. Changing lead geometry needs a matching new FEM solve.

`sample_fem.py` runs in the existing FEM environment, loads the large mesh once
and writes a small `sampled_field.npz`. NEURON only loads those sampled values.
All requested compartments must be located in a tetrahedron; **no colour clamp,
nearest-node fill, or silent extrapolation** is used. Units are explicit:
coordinates mm → µm for axon geometry, `V/A × mA × 10^-3 × 10^3 = mV`.

The supplied polyline is a strictly sampled **white-matter engineering fixture**
for the currently cached lead. It is not an anatomically justified Aδ or C
trajectory and must not be used as evidence of their recruitment in that tract.

## Interpretation and remaining work

The detector requires outward 0-mV crossings at two separated recording sites;
simultaneous local excursions do not qualify. Inspect the saved nine-site traces.
The doubling/bisection search assumes a monotonic local excitation bracket;
it is not a general search for activation/block windows at extreme amplitudes.
Recording locations and sealed-end effects must be checked for actual paths.

[docs/neuron_plan.md](../../docs/neuron_plan.md) separates the near-term poster
scope from later grant work. Remote SSH/SLURM remains a TODO. The simulator and
field adapter are command-line modules; a NEURON FreeCAD macro has not been added.
