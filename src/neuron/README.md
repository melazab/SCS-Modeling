# `src/neuron/` — cable modelling

Takes the extracellular potential field from the Ansys solve and works out which
nerve fibres it actually activates, producing recruitment curves for dorsal
versus ventral lead placement.

Read [`docs/neuroanatomy_primer.md`](../../docs/neuroanatomy_primer.md) first if
the neuroscience is unfamiliar — it explains what is being stimulated and why,
and it is what makes the modelling choices below make sense. The staged plan is
in [`docs/neuron_plan.md`](../../docs/neuron_plan.md).

## Install

    python3 -m venv src/neuron/.venv
    src/neuron/.venv/bin/pip install -r src/neuron/requirements.txt
    src/neuron/.venv/bin/python src/neuron/smoke_test.py

Its own venv, on purpose — NEURON ships a compiled extension plus its own
binaries and wants to own its environment. The other two directories use the
system interpreter; this one does not.

## Status

| chunk | what | state |
|---|---|---|
| 0 | environment + smoke test | **done** |
| 1 | MRG axon, validated by conduction velocity | next |
| 2 | analytic point source, strength–duration | |
| 3 | Ansys → NEURON field export contract | |
| 4 | axon trajectories from the anatomy | |
| 5 | coupling, activating-function check | |
| 6 | single-fibre threshold bisection | |
| 7 | population sweep (SLURM array) | |
| 8 | dorsal vs ventral comparison | |

Chunks 0–2 need no FEM at all and run on a laptop in minutes. That is
deliberate: they build and validate the cable machinery independently of whether
the big solve is behaving, so a later disagreement can be localised to one half
of the pipeline.

## `smoke_test.py`

Checks two things that fail independently:

1. The `neuron` module imports and can integrate a cable equation — a
   Hodgkin–Huxley soma is injected with current and must spike. HH is not a model
   of anything in this study (it is unmyelinated squid at 6.3 °C); it is just the
   cheapest end-to-end exercise of section, mechanism, stimulus and integrator.
2. **`nrnivmodl` compiles a `.mod` file into a loadable mechanism.** This needs a
   working C compiler, and it is the one that matters — the MRG axon in chunk 1
   ships its channel mechanisms as `.mod` files. A NEURON that imports but cannot
   compile is useless here, and without this check the failure would not surface
   until much later.

One trap, since it looks exactly like a compile failure: NEURON **auto-loads**
`x86_64/libnrnmech.so` from the working directory at import. Calling
`h.nrn_load_dll()` on it as well raises `The user defined name already exists`,
which is in fact proof the compile succeeded.

## The fibre model, and why not Hodgkin–Huxley

Chunk 1 builds the **MRG** axon (McIntyre, Richardson & Grill 2002, ModelDB
#3810): a double-cable myelinated fibre with explicit nodes of Ranvier, paranodal
and internodal segments, and a periaxonal space. That structure is not decoration
— extracellular stimulation drives a fibre through the *second spatial derivative*
of the potential along it, and where the nodes fall relative to the field is what
sets threshold. HH has no nodes and would answer the wrong question.

It is not a pip package: download the `.mod` files from ModelDB and compile them
with `nrnivmodl`.
