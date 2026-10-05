# Representative Aβ, Aδ and C axons

These models establish a reproducible three-fibre pipeline for preliminary
poster work. They are **isolated axon adaptations**, not whole sensory neurons,
and are not interchangeable representations of all sensory subtypes.

## Source selection

| Label | Implementation | Initial geometry | Temperature used here |
|---|---|---|---|
| Aβ reference | MRG 2002 large myelinated mammalian axon | 5.7 µm outer fibre diameter; exact published geometry table entry | 37 °C |
| Aδ reference | Dorsal-root portion of Hao et al. 2023 / Jaffe DRGsims | 3 µm nodes and internodes, 1.5 µm nodes, 150 µm internodes | 37 °C |
| C reference | Peripheral axon of Sundt et al. 2015, Figure 5d | 0.8 µm unmyelinated axon | 37 °C (upstream example uses 35 °C; CLI supports that too) |

Sources:

- McIntyre, Richardson & Grill (2002), *Modeling the excitability of mammalian
  nerve fibers: influence of afterpotentials on the recovery cycle.*
  [ModelDB 3810](https://modeldb.science/3810),
  [source](https://github.com/ModelDBRepository/3810).
- Hao et al. (2023), *Dorsal root ganglia control nociceptive input to the central
  nervous system.* [Paper](https://doi.org/10.1371/journal.pbio.3001958),
  [author's source](https://github.com/dbjaffe67/DRGsims).
- Sundt, Gamper & Jaffe (2015), *Spike propagation through the dorsal root ganglia
  in an unmyelinated sensory neuron: a modeling study.*
  [Paper](https://doi.org/10.1152/jn.00226.2015),
  [ModelDB 187473](https://modeldb.science/187473).

`vendor/sources.json` records exact upstream commits, paths and SHA-256 hashes.
Vendored channel files and reference HOC files are unmodified. The build checks
their hashes before compilation. Attribution and source README files are kept
with them; this project does not relicense those upstream files.

## Adaptations and limits

- MRG has active nodes, MYSA/FLUT/STIN compartments, and a periaxonal cable.
  The source's membrane properties are retained. Its original one-segment-per-
  section discretization showed a roughly 20% threshold change when subdivided
  in the analytic extracellular test. This implementation therefore uses three
  segments per section and checks against five; it is **not a byte-equivalent
  reproduction of the original discretized model**. Published field/threshold
  replication remains a separate task.
- The Aδ implementation uses the source's central axon parameters: nodal TTX-
  sensitive sodium current, passive internodes, and resting-current leak balance.
  It omits the soma, stem, T-junction and GABA conductance. No potassium channel
  was added to the nodes: the source does not put one there. This is a simplified
  representative Aδ model, not a validated human nociceptor subtype model.
- The C implementation retains the source axon's shifted sodium and delayed-
  rectifier potassium currents, passive leak, and resting-current balance.
  DRG-local M-current, soma and T-junction are omitted. This model's recovery
  cycle and response to sustained high-frequency trains are not validated here.
- No diameter distributions, branching, population recruitment, synapses or
  analgesia model are implemented. Aβ/Aδ/C labels identify the intended initial
  representatives, not a claim that diameter alone determines sensory function.
- All three have sealed ends. Test length and recording-site sensitivity before
  interpreting anatomical thresholds; simultaneous voltage excursions are not
  counted as propagation.

## Validation levels

The automated checks establish mechanism compilation, resting stability,
intracellular propagation, zero/uniform extracellular-field controls, strict
units/contact mapping, and numerical threshold refinement. They do not establish
physiological validation by themselves. The next scientific checks are published
waveform/CV/recovery comparisons, strength–duration curves, and anatomically
reviewed trajectories for each fibre class.

The supplied FEM trajectory is deliberately labelled an **engineering fixture**
in white matter. Applying all three models there tests the software coupling;
it does not assert that Aδ and C dorsal-root afferents follow dorsal-column paths.
