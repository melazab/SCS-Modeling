# NEURON plan — poster milestone and research extension

## Scope update — September 30, 2026

The immediate poster target is now a matched dorsal/ventral FEM field comparison.
NEURON and downstream neuronal-population modeling are deferred from the minimum
poster deliverable. The existing three-axon software milestone remains available;
the longer-term representatives are Aβ, Aδ and C fibres, not three anatomical
populations of the same MRG model. See [poster_outline.md](poster_outline.md)
and [the post-poster FEM roadmap](../fem/TODO.md). The sequence below describes
subsequent axon work, not completed poster results.

## Current software milestone

Implemented under `src/neuron/`: published-source axon adaptations, compiled
mechanisms, propagation checks, analytic extracellular stimulation, balanced
biphasic pulses, propagated-spike threshold search, numerical refinement checks,
and strict sampling of cached Elmer basis fields at every compartment.
See [README](../src/neuron/README.md) and
[model provenance and limits](../src/neuron/MODEL_SOURCES.md).

The earlier plan assumed Ansys, eight contacts and a grounded outer boundary.
The production implementation now uses Elmer, the actual contact count, a common
distributed return and a reference potential. Balanced contact combinations
cancel the shared return. Existing fields are in V/A at coordinates in mm.
A waveform combines those spatial fields with time-varying contact currents;
NEURON receives mV at its compartments. The coupling is one-way.

## Subsequent axon-modeling sequence

1. **Establish three functioning representative axons.** Reproduce source
   membrane behaviour and propagation, check rest stability, numerical
   convergence and extracellular coupling. Distinguish software checks from
   experimental/physiological validation. The Aδ source is a simplified central
   axon; the C source is unmyelinated and not HH squid. Model choice/adaptation
   must accompany preliminary figures.
2. **Validate stimulation responses.** Compare published AP/CV and recovery
   behaviour, generate strength–duration curves, check polarity, and test axon
   length and node-position sensitivity. Expected trends are checks to
   investigate, not rules enforced by retuning channels.
3. **Review anatomical trajectories.** Define dorsal-root paths for the three
   afferent types and appropriate central continuations. Do not put Aδ and C
   fibres along dorsal-column paths merely to make all three use one line.
   Render paths in FreeCAD for review; assert coordinate transforms, laterality,
   tissue membership and full FEM coverage. The current white-matter polyline
   is explicitly a software fixture only.
4. **Small matched dorsal/ventral comparison.** Hold fibre parameters, paths,
   pulse definition and contact montage conventions fixed across placements.
   Report thresholds or no activation at tested amplitudes, propagation traces,
   and the limits of three representatives. These are not population recruitment
   curves. Reserve time for figures and repeatability rather than a large sweep.

## FEM dependencies before interpreting anatomical thresholds

- White-matter anisotropy and root Middle-layer conductivity remain unresolved.
- Check field/threshold convergence along the actual axon paths, including
  piecewise-linear interpolation noise and boundary/truncation sensitivity.
- Changing radial contact projection changes the contact outer surface;
  regenerate the mesh and basis fields and resolve the resulting surface step.
- Earlier field agreement with RADO examples is not validation of the current
  geometry, material assumptions or neuronal threshold predictions.

## Broader R01 work after the first milestone

Population diameters and trajectories; dorsal-column and spinothalamic targets;
ventral-root motor recruitment as a competing threshold; tonic/burst/10-kHz
trains with entrainment and conduction failure; root/DRG branching and channel
heterogeneity; downstream circuit models for indirect effects; and remote
SSH/SLURM sweeps. Reconcile T6–T11 study coverage with available anatomy.
Direct axon firing alone does not establish analgesia or circuit modulation.
