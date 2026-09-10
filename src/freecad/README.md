# `src/freecad/` — geometry side

Everything that reads or writes the STL set and the `NBF_RADO-SCS.FCStd`
document: labelling bodies, measuring the epidural space, building parametric
SCS leads, and colouring the model to match Ansys Engineering Data.

## Install

    pip install -r requirements.txt

and a FreeCAD ≥ 1.0 install providing `freecadcmd` on PATH — see the notes at the
top of `requirements.txt`, since FreeCAD itself cannot come from pip. Developed
against FreeCAD 26.3.0 / Python 3.14.4.

## Two kinds of script

| script | interpreter | what it does |
|---|---|---|
| `apply_labels.py` | `freecadcmd` | gives all 245 bodies human-readable `Label`s from `body_aliases.yaml` |
| `apply_colors.py` | `freecadcmd` | colours bodies by tissue, using `../ansys/tissue_map.yaml` |
| `measure_corridor.py` | `freecadcmd` | ray-casts the epidural mesh, writes `epidural_corridor.json` |
| `make_scs_lead.py` | `freecadcmd` | builds an *n*-contact lead swept along the measured canal centreline |
| `build_lead_config.py` | `freecadcmd` (`--list`/`--show` work anywhere) | builds and **validates** any named configuration from `lead_configs.yaml` |
| `make_tissue_groups.py` | running FreeCAD GUI | puts every body in a group named for its tissue |
| `check_laterality.py` | `python3` | re-derives left/right from geometry and validates `body_aliases.yaml` |
| `find_floating_bodies.py` | `python3` | proximity graph over the STL surfaces; finds electrically isolated islands |

### Running the `freecadcmd` ones

`freecadcmd` owns the command line — it consumes every flag itself, and anything
after its `--pass` makes it skip the script entirely. So options go through an
environment variable instead:

    freecadcmd src/freecad/measure_corridor.py
    MAKE_LEAD_ARGS="--side ventral --z-center 110" freecadcmd src/freecad/make_scs_lead.py
    APPLY_LABELS_ARGS="--dry-run" freecadcmd src/freecad/apply_labels.py
    BUILD_LEAD_CONFIG_ARGS="--all --validate" freecadcmd src/freecad/build_lead_config.py

Two more `freecadcmd` behaviours these scripts are written around: it *imports*
the script rather than running it as `__main__` (so a `__main__` guard never
fires, and the work happens at import time), and it swallows stdout (so every
script also writes a report file next to its output). It also segfaults on exit
*after* the document is saved, which is harmless.

## Why the leads have to follow a curve

The T8–T10 spine here is kyphotic, so the epidural space migrates posteriorly as
you go rostral — about 7 mm over the 43 mm length of an 8-contact lead. A lead at
fixed *y* would walk straight out of the space and through the dura. So
`measure_corridor.py` ray-casts the epidural mesh to fit a centreline *y(z)*
(cubic, 0.018 mm RMS over 98 mm) and `make_scs_lead.py` sweeps the lead along it.

Ray casting, not surface arithmetic: RADO's compartments are **hollow shells that
tile space**, each with the inner ones carved out. Measuring a channel as
`epidural.ymax − dura.ymax` is the difference of two *outer* surfaces and is not
the thickness of anything — it once produced a confident and completely wrong
"there is no dorsal epidural space" conclusion. Casting a ray along +Y and
reading consecutive entry/exit pairs gives the real numbers:

    dorsal epidural fat   2.27 – 2.35 mm
    ventral epidural fat  1.67 – 1.74 mm

so a 1.3 mm clinical lead fits on either side with no carving of the anatomy.

## Lead configurations

A lead is described clinically — “8 contacts, 3 mm long, 1 mm apart, dorsal, at
T10, 1 mm left of midline” — and `lead_configs.yaml` is that description, one
block per named configuration, inheriting from a `defaults:` block.
`build_lead_config.py` turns a name into geometry and, more usefully, into a
verdict:

    python3 src/freecad/build_lead_config.py --list
    BUILD_LEAD_CONFIG_ARGS="--name khadka_clinical_T10 --validate" \
        freecadcmd src/freecad/build_lead_config.py
    BUILD_LEAD_CONFIG_ARGS="--name dorsal_T11_left1 --export" \
        freecadcmd src/freecad/build_lead_config.py

Any field can be overridden on top of a named entry (`--level T10 --x-offset -2`),
which is how the GUI drives it. `src/freecad/lead_designer.FCMacro` is the same
thing as a dock panel: spinboxes, a green/red verdict with the measurements, a
disposable preview in the open document, and an Export button. Install it like
`tissue_visibility.FCMacro`.

### Three things it checks, and why there is no boolean subtraction

The intuition is that a lead must be *carved* out of the tissue around it. It
must not, for two reasons. RADO's compartments are hollow shells that tile
space, and the epidural fat is a real void the lead occupies — there is nothing
to subtract. And `Mesh.difference()` is a **complete no-op in this FreeCAD
build**, verified: it returns cleanly and changes nothing, so anything built on
it silently produces an unmodified model.

What replaces carving is a fit check, in three measurements:

| check | what it measures |
|---|---|
| fit | channel thickness ray cast **at that side, offset and z span**, against the diameter |
| containment | every sampled vertex of the built lead: inside the epidural fat, outside the dura |
| coverage | whether the lead runs past the end of the measured centreline |

A configuration that fails is reported in clinical terms and **not exported** —
“the lead is 2.00 mm across, but the ventral epidural fat is only 1.68 mm thick
at its narrowest; it would press through the dura.” `ventral_z110_2mm_FAILS` in
the catalogue is kept as a worked example of exactly that.

Lateral offsets get their **own** centreline. The stored polynomial in
`epidural_corridor.json` is measured at the midline, and the channel centre
migrates about 1.1 mm of *y* by 4 mm off-midline, so `build_lead_config.py`
re-ray-casts and re-fits at the requested *x*. At the midline it uses the stored
fit and prints how far the two disagree, so they cannot drift apart unnoticed.

### Vertebral levels

`levels:` in the YAML maps `T10` / `T11` / `T12` (and the four disc levels) to a
*z*, derived from the geometry: RADO names its discs for the levels they
separate, the four come out in the order “+Z is rostral” demands, and a
vertebral body is the gap between two of them. Cross-checked against the
anterior mass of each vertebra STL — the two derivations agree to 0.6 mm.

**Caveat, and it is not a small one:** those disc names imply the three bodies
are T10–T12, while the paper says T9–T11 and RADO's filenames say “T8-10”. All
three disagree and nothing in the geometry settles it. The *z* coordinates are
solid; the level *names* are RADO's, taken at face value.

## Generated output, and why there is no `.FCStd` per configuration

`generated_leads/<tag>/` holds one STL per contact plus one for the insulator,
already positioned in the model's own coordinates — import them without any
transform. `dorsal_z110_8c/` and `ventral_z110_8c/` are the current pair; both
verified at 100% of sampled vertices inside epidural fat and 0% in dura, and
both reproduced byte-for-byte by `build_lead_config.py --export`.

The document is 25 MB, so ten configurations must not mean ten documents, and
every save of one is another full copy in git history. The durable artefact is
therefore the YAML entry plus its STL directory — a dozen numbers and a few
hundred kB. Bodies in a FreeCAD document are a **disposable preview**: the macro
(or `build_lead_config.preview()`) drops the lead into the active document
inside a group named `SCS_LeadPreview`, replacing any previous preview, and
“Remove preview” puts the tree back. Explore, export the one you want, close the
document without saving.

`NBF_RADO-SCS_dorsal.FCStd` and `_ventral.FCStd` are the exception, not the
pattern: they are the two arms of the dorsal-vs-ventral study and are what Ansys
is handed. They were built by `make_lead_variants.py` and they are the only two
that should exist.

**Nothing in this pipeline saves a model document headlessly.** A `.FCStd`
written with no GUI running has no `GuiDocument.xml` and no `ShapeAppearance`
blobs at all — every tissue colour gone, and the document still opens looking
like a fresh grey import. `make_tissue_groups.py` documents the measurement;
`build_lead_config.preview()` refuses outright to touch a saved document when
`FreeCAD.GuiUp` is false.
