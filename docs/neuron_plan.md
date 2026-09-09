# NEURON stage — plan

Goal: turn the FEM potential field into **recruitment curves** — what fraction of
each fibre population is activated, at what amplitude — so ventral and dorsal
epidural stimulation can be compared at several lower-thoracic levels.

## The shape of the whole thing

The FEM and the cable model are coupled in one direction only, and that is what
makes this tractable:

    Ansys                          NEURON
    -----                          ------
    solve conduction   ->  phi(x)  ->  sample phi at each node of Ranvier
    (slow, hours)                      run the axon (fast, milliseconds)
                                       bisect amplitude -> threshold
                                       repeat over a fibre population

Quasi-static conduction is **linear**, so phi scales with amplitude and
superposes across contacts. Nothing in NEURON feeds back into the FEM.

### The one decision that has to be made early

**Solve once per contact, not once per stimulation pattern.**

For each of the 8 contacts, run one solve with 1 A into that contact and the
outer boundary grounded, giving `phi_i(x)` per unit current. Then *any*
configuration — bipolar, tripolar, guarded cathode, any amplitude, any current
steering ratio — is just

    phi(x) = sum_i  I_i * phi_i(x),      with sum_i I_i = 0

computed in milliseconds in NEURON. Eight solves buy you the entire parameter
space. Solving per configuration instead would mean a fresh multi-hour `smp` job
for every point on every recruitment curve, which is simply not affordable.

This has to be settled before the production sweep starts, because it changes
what the decks do. Everything else below can be reordered.

## Chunks

Chunks 0–2 need **no FEM at all** and run on the laptop in minutes. Start there —
they are independent of whether the big Ansys solve is behaving, and they build
the thing that will otherwise be impossible to debug later.

---

### Chunk 0 — environment (half a day)

`pip install neuron` (wheels exist for 8.2+), confirm `nrnivmodl` compiles a mod
file, run one Hodgkin–Huxley soma and see a spike.

Deliverable: `src/neuron/` with `requirements.txt`, `README.md`, one smoke test.

### Chunk 1 — the axon, in isolation (1–2 days)

Build the **MRG** double-cable myelinated axon (McIntyre, Richardson & Grill
2002; ModelDB 3810). It is the standard fibre model for SCS work — 11 segments
per internode, explicit myelin and periaxonal space, so it responds to
extracellular fields the way a real myelinated fibre does. Diameters 5.7–16 µm.

**Validation, and this is the point of the chunk:** measure conduction velocity
vs diameter and check it against the published table (roughly CV ≈ 4.6·D m/s, so
~26 m/s at 5.7 µm to ~75 m/s at 16 µm). If CV is right, the axon is right.

Deliverable: `mrg_axon.py` (build an axon of N nodes at a given diameter),
`test_conduction_velocity.py`.

### Chunk 2 — analytic point source (1–2 days)

Stimulate that axon with a monopole in an infinite homogeneous medium, where the
answer is known in closed form: `phi = I / (4·pi·sigma·r)`.

Reproduce three textbook results: the strength–duration curve (chronaxie ~100–200
µs for large myelinated fibres), threshold rising steeply with electrode–fibre
distance (roughly as distance²), and cathodic threshold below anodic.

This validates the **extracellular coupling machinery** independently of the FEM.
Do it now and later, when the full pipeline gives a strange number, you will know
which half to suspect. Skip it and you will not.

Deliverable: `point_source.py` and a strength–duration figure.

### Chunk 3 — the Ansys → NEURON contract (2–3 days)

Freeze the interface. This is the chunk that constrains the FEM work, so do it
before the production solves.

- APDL `/POST1` snippet that writes `phi` per unit current on a structured grid
  covering cord + roots, ~0.1–0.25 mm spacing, plus a variant that writes `phi`
  at an arbitrary list of points (for exact node coordinates).
- A documented `.npz` schema: grid origin, spacing, shape, `phi_i` for each
  contact, units, and the contact geometry it came from.
- Sanity checks baked in: `phi` must fall as ~1/r far from the contact, and
  current in must equal current out.

Watch the units at this boundary — the decks are `/units,uMKS` (µm, pA,
Tohm·µm) and NEURON wants mV, µm, ms. Convert once, in `field.py`, and assert it.

Deliverable: `src/ansys/export_field.dat`, `src/neuron/field.py` (load +
trilinear interpolation).

### Chunk 4 — axon trajectories from the anatomy (2–3 days)

Where do the fibres actually go? Generate 3D polylines from the existing
geometry, reusing the ray-casting already in `src/freecad/measure_corridor.py`:

- **Dorsal column fibres** — longitudinal in dorsal white matter, sampled over a
  range of dorsoventral depths and mediolateral offsets. These are the intended
  target of dorsal SCS.
- **Dorsal root fibres** — following the root bodies from the DRG into the dorsal
  horn. These typically have *lower* thresholds than dorsal columns, which is why
  they matter clinically.
- **Ventral root / motor fibres** — the ones ventral stimulation might recruit,
  and the reason the ventral arm of the study is interesting at all.

Deliverable: `src/freecad/make_axon_paths.py` writing polylines to JSON, plus a
render of the paths inside the cord so they can be eyeballed before anything is
run on them.

### Chunk 5 — glue (1–2 days)

Place MRG nodes along a trajectory by arc length, sample `phi` at each, and apply
it via the standard `xtra` + `extracellular` pattern.

**Go/no-go check:** plot `phi` along one axon and its second spatial difference —
the activating function. It should be smooth, and should peak under the cathode.
If it is noisy, the FEM grid is too coarse; if it peaks somewhere odd, the
trajectory or the transform is wrong. Do not proceed past a bad-looking plot.

Deliverable: `couple.py` and that figure.

### Chunk 6 — single-fibre threshold (1 day)

Bisection on amplitude, AP detection at a distal node, with a biphasic pulse of
the clinical pulse width. Returns a threshold in mA for a given (fibre,
trajectory, contact configuration, waveform).

Deliverable: `threshold.py`.

### Chunk 7 — population sweep (2–3 days)

Sweep fibre diameter × trajectory × dorsal/ventral × rostro-caudal level. Every
run is independent, so this is a SLURM **array** job — and being NEURON rather
than Ansys it is cheap, so it does not collide with the 24-CPU `tlv` cap the way
the FEM does.

Deliverable: recruitment curves — % of population activated vs amplitude, one
curve per fibre class per lead placement.

### Chunk 8 — the actual comparison

Ventral vs dorsal, per level, per fibre class. The figures for the paper.

---

## Validation checkpoint before believing anything

With the dorsal lead in place, dorsal **root** fibres should recruit at lower
amplitudes than dorsal **column** fibres, and absolute thresholds should land in
the range published for percutaneous SCS. Reproducing that ordering is the
strongest available check that the whole chain is sound. If it inverts, suspect
Chunk 4 (trajectories) or Chunk 3 (units) before suspecting the biophysics.

## Risks

- **Grid resolution.** The activating function is a *second difference* of `phi`,
  so it amplifies interpolation error. Chunk 5's plot is the early warning.
- **Chunk 3 before the production solves.** Getting this backwards means
  re-running multi-hour jobs to export something forgotten.
- **The white matter is currently isotropic.** Conductivity anisotropy (0.1432
  transverse / 0.6 longitudinal S/m) has a direct effect on longitudinal fibre
  recruitment. Fine to defer, but it must be revisited before any claim about
  absolute thresholds — see `src/ansys/README.md`.
- **The root "Middle" layer question** is still open (CSF vs nerve, a 12×
  conductivity difference right next to the electrodes). It should be settled
  before Chunk 7, not after.
