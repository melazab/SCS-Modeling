# `src/freecad/` — geometry side

Everything that reads or writes the STL set and the `NBF_RADO-SCS.FCStd`
document: labelling bodies, measuring the epidural space and the neuroforamina,
building parametric SCS leads, and colouring the model to match Ansys
Engineering Data.

**There is one model document.** `NBF_RADO-SCS.FCStd` holds anatomy and nothing
else. A lead is not baked into a copy of it — it is an entry in
`lead_configs.yaml` plus a directory of STLs, previewed into the open document
on demand and thrown away. The former `NBF_RADO-SCS_dorsal.FCStd` and
`_ventral.FCStd`, and the `make_lead_variants.py` that built them, are gone.

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
| `measure_foramen.py` | `freecadcmd` | ray-casts all eight neuroforamina, writes `foraminal_corridors.json` |
| `make_scs_lead.py` | `freecadcmd` | builds an *n*-contact lead swept along a measured centreline |
| `build_lead_config.py` | `freecadcmd` (`--list`/`--show` work anywhere) | builds and **validates** any named configuration from `lead_configs.yaml` |
| `make_tissue_groups.py` | running FreeCAD GUI | puts every body in a group named for its tissue |
| `drop_rado_lead.py` | running FreeCAD GUI | removes RADO's own 4-contact DRG lead from the anatomy document |
| `check_laterality.py` | `python3` | re-derives left/right from geometry and validates `body_aliases.yaml` |
| `find_floating_bodies.py` | `python3` | proximity graph over the STL surfaces; finds electrically isolated islands |

### Running the `freecadcmd` ones

`freecadcmd` owns the command line — it consumes every flag itself, and anything
after its `--pass` makes it skip the script entirely. So options go through an
environment variable instead:

    freecadcmd src/freecad/measure_corridor.py
    freecadcmd src/freecad/measure_foramen.py
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

## And why a DRG lead follows a different curve entirely

A percutaneous **DRG** lead does not go down the canal at all. It leaves the
epidural space at one level and runs *laterally*, out through the neuroforamen,
so its contact array ends up beside one dorsal root ganglion. Its long axis is
*x*, and both of the other coordinates have to follow a curve — so
`measure_corridor.py`'s *y(z)* at fixed *x* is the wrong shape of answer.

`measure_foramen.py` measures the right one. RADO models the tissue filling each
foramen as one **solid** block (`additional_menginges-2..9`, four a side, one
per foramen, each wrapped around a root and its ganglion). Unlike the canal
compartments they are not hollow shells, which is physically right: a DRG lead is
threaded *through* foraminal fat, displacing it, not dropped into a void. The
block's centre at each lateral station *u* is found by **alternating** ray casts,
seeded on the ganglion and iterated to a fixed point — cast +Z to get *z* and the
z-extent, cast +Y at that *z* to get *y*, repeat — and *y(u)*, *z(u)* are fitted
as quartics over the run where the block is a real channel (both extents ≥ 5 mm).

All eight measure cleanly: about 20 mm of usable channel each, at least
5.2 × 6.1 mm of cross-section, quartic residuals ~0.13 mm RMS in *y* and
0.07–0.15 mm in *z*, worst case ~0.5 mm. **That is an order of magnitude looser
than the canal fit** (0.018 mm RMS) and the reason is real: the foramen is a
short irregular channel between two pedicles, not a smooth canal. Read ~0.5 mm as
the honest positional uncertainty of a generated DRG lead.

### Is that trajectory right? RADO's own lead says mostly yes

RADO ships a 4-contact lead 11–25 mm left of midline at *z* ≈ 93–97, on a curved
insulator: a DRG lead, on the left level-3 ganglion, placed **by hand in
SolidWorks** and derived from nothing in this pipeline. It is the only
independent check available, and `--compare-rado` measures the agreement rather
than asserting it. Sampling the measured corridor at RADO's own four contact
positions:

| at *u* | across-corridor gap |
|---|---|
| 11.16 mm | 4.27 mm |
| 15.65 mm | 0.58 mm |
| 20.30 mm | 1.39 mm |
| 25.40 mm | 1.32 mm |

**Over the ganglion and distal to it the corridor reproduces RADO's lead to about
1.1 mm; medial of the ganglion it does not, and by a lot.** The reason is
visible in the measurement: medial of the ganglion the block flares vertically to
meet the canal — 14 mm tall in *z* at 11 mm out — and its geometric centre climbs
3–4 mm above the nerve root, while RADO's lead stays down on the root. So the
corridor is a good description of the foramen proper and a poor one of the
lateral recess. That is why `lateral_offset` defaults to 1.0 mm rather than 0:
it keeps every contact out where the corridor is trustworthy.

Two more things the comparison settled, both of which changed the checks:

* RADO's own lead is **100% inside the foraminal tissue, 0% in the ganglion,
  0% in the dura** — so those are the right checks, and they are what a DRG lead
  is validated against.
* RADO's own lead **does** clip three nerve-root bodies (7, 8 and 56 sampled
  points), because the foraminal block is not carved out around the roots. So
  grazing a root is a WARN here, not a FAIL. Only the ganglion and the thecal
  sac are hard failures.

### The standoff is solved for, not written down

The centre of the foraminal corridor runs along the **top of the ganglion** — on
L3 it passes within about 0.1 mm of the surface, and RADO's hand-placed lead sits
in the same fraction of a millimetre. That is clinically correct (a DRG lead is
supposed to lie against its ganglion) but it means a lead dropped exactly on the
corridor centre is a coin toss between grazing it and cutting into it.

So the clinical parameter is `ganglion_clearance` — “hold the lead this far off
the ganglion” — and the rostral shift needed to achieve it is **computed per
lead**. It has to be: 0.25 mm of clearance needs +0.77 mm on L3 and +1.62 mm on
R3, and either number written into the catalogue would be wrong for the other
side. Rostrally rather than dorsally because that is where the room is (the block
is 13–14 mm tall in *z* against 8–10 mm in *y*, with the ganglion low in it) and
because that is where a percutaneous DRG lead goes in theatre: into the superior
aspect of the foramen, under the pedicle, over the ganglion.

## Lead configurations

A lead is described clinically — “8 contacts, 3 mm long, 1 mm apart, dorsal, at
T10, 1 mm left of midline”, or “four contacts out to the left third ganglion” —
and `lead_configs.yaml` is that description. `build_lead_config.py` turns a name
into geometry and, more usefully, into a verdict:

    python3 src/freecad/build_lead_config.py --list
    BUILD_LEAD_CONFIG_ARGS="--name dorsal_pair_T11 --validate" \
        freecadcmd src/freecad/build_lead_config.py
    BUILD_LEAD_CONFIG_ARGS="--name rado_drg_L3_match --validate --compare-rado" \
        freecadcmd src/freecad/build_lead_config.py
    BUILD_LEAD_CONFIG_ARGS="--name dorsal_T11_left1 --export" \
        freecadcmd src/freecad/build_lead_config.py

Any field can be overridden on top of a named entry (`--level T10 --x-offset -2`,
`--target R2`), which is how the GUI drives it. `src/freecad/lead_designer.FCMacro`
is the same thing as a dock panel: a tab of spinboxes per lead, a green/red
verdict with the measurements, a disposable preview in the open document, and an
Export button. Install it like `tissue_visibility.FCMacro`.

### Three lead types

| `type:` | runs | positioned by | swept along |
|---|---|---|---|
| `dorsal` | rostro-caudally, dorsal canal | `level`/`z_center`, `x_offset` | `epidural_corridor.json` |
| `ventral` | rostro-caudally, ventral canal | `level`/`z_center`, `x_offset` | `epidural_corridor.json` |
| `drg` | laterally, out through a foramen | `target` (`L1`..`L4`, `R1`..`R4`), `lateral_offset`, `ganglion_clearance` | `foraminal_corridors.json` |

They are not three settings of one thing, so their parameters do not overlap and
`resolve_lead()` **rejects** rather than ignores one that means nothing for the
type — `z_center` on a DRG lead, or `target` on a dorsal one, is an error, on the
same principle that a typo'd `diamter: 2.0` is.

`side: dorsal|ventral` is still accepted as a spelling of `type:`, for older
entries and the `--side` flag. It cannot spell `drg`: a DRG lead is not a side of
the canal.

**Ganglion names are the model's own, not lumbar levels.** `L3` means “left,
third from the top”, counted rostral to caudal from the ganglion cores' *z* —
the same derivation `body_aliases.yaml` uses for its labels, and it reproduces
that file's `level_map` exactly. In this T8–T10 model `L3` is a *thoracic*
ganglion at about *z* = 93.

### A configuration is a set of leads

Clinicians do not implant one lead and stop. Two dorsal leads straddling the
midline is routine; so is a dorsal plus a DRG lead, or one DRG lead per side. So
an entry is a **set**, written either way:

```yaml
dorsal_z110_8c:              # one lead: describe it inline, as before
  z_center: 110.43203353881836

dorsal_pair_T11:             # several: a `leads:` list
  level: T11                 #   keys out here apply to EVERY lead in the set
  leads:
    - x_offset:  1.5
      contacts:  8
    - x_offset: -1.5
```

`max_leads_per_type` at the top of the file caps how many of any **one** type a
set may hold; it is 2 by default, so the largest legal set is two dorsal + two
ventral + two DRG. The Add buttons in the macro stop there too.

Layering for one lead of a set, later wins:

    defaults  ->  epidural_defaults | drg_defaults (by type)  ->  the
    configuration's shared keys  ->  that lead's own keys  ->  CLI/GUI overrides

An entry with no `leads:` list is a set of one whose single lead is the entry
itself, which is why every pre-existing entry still resolves to exactly the
geometry it always did.

### Three things it checks, and why there is no boolean subtraction

The intuition is that a lead must be *carved* out of the tissue around it. It
must not, for two reasons. RADO's compartments are hollow shells that tile
space, and the epidural fat is a real void the lead occupies — there is nothing
to subtract. And `Mesh.difference()` is a **complete no-op in this FreeCAD
build**, verified: it returns cleanly and changes nothing, so anything built on
it silently produces an unmodified model.

What replaces carving is a fit check. For a **dorsal or ventral** lead, three
measurements:

| check | what it measures |
|---|---|
| fit | channel thickness ray cast **at that side, offset and z span**, against the diameter |
| containment | every sampled vertex of the built lead: inside the epidural fat, outside the dura |
| coverage | whether the lead runs past the end of the measured centreline |

A **DRG** lead asks a different question, because the foramen is not a void:

| check | what it measures |
|---|---|
| fit | clearance from the lead axis to the wall of the foraminal block, in *y* and *z*, at every station |
| containment | every sampled vertex inside that block |
| clearance | every sampled vertex **outside** the target ganglion and outside the thecal sac (FAIL); overlap with a nerve root is reported but only WARNs — RADO's own lead does it too |
| coverage | whether the array runs off the fitted corridor (WARN — the containment check is what decides) |

A configuration that fails is reported in clinical terms and **not exported** —
“the lead is 2.00 mm across, but the ventral epidural fat is only 1.68 mm thick
at its narrowest; it would press through the dura.” `ventral_z110_2mm_FAILS` in
the catalogue is kept as a worked example of exactly that.

### And one check no single lead can make

Two leads can each pass everything above and still be un-implantable, because
nothing in a per-lead check looks at the other lead. So the **set** is checked
too: every pair, by exact shape-to-shape distance, and when that is zero, by the
volume they share and how much of one lead's sampled surface is inside the other.

    dorsal_pair_T11         ±1.5 mm off midline → 1.70 mm apart      PASS
    dorsal_pair_tight_FAILS ±0.5 mm off midline → 0.00 mm, INTERSECT FAIL

`dorsal_pair_tight_FAILS` is kept for the same reason `ventral_z110_2mm_FAILS`
is: both its leads fit the channel perfectly well on their own.

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
transform. `dorsal_z110_8c/` and `ventral_z110_8c/` are what the two former study
documents contained; both verified at 100% of sampled vertices inside epidural
fat and 0% in dura, and both **reproduced byte-for-byte** by
`build_lead_config.py --export`. That is this module's regression test: any
change to the resolve/sweep path has to leave all eighteen files identical.

The layout is deliberately not uniform:

    one lead    generated_leads/<tag>/SCS Lead Electrode 1.stl, ...
    several     generated_leads/<tag>/lead1_dorsal/SCS Lead Electrode 1.stl, ...
                generated_leads/<tag>/lead2_drg_L3/...
                generated_leads/<tag>/leads.txt        (which lead is which)

A single-lead configuration keeps the flat layout it has always had, because
those two directories are committed at those exact paths and moving them would
make the byte-for-byte test untestable. A set gets one subdirectory per lead,
because the **filenames** cannot carry the distinction: they are RADO's own, and
`../ansys/tissue_map.yaml` matches the object `Name`s they sanitise to, so
renaming them costs the contacts their silver and the insulator its material.
The directory carries it instead — and note that importing two leads into one
FreeCAD document therefore gives the second one's bodies uniquified `Name`s
(`SCS_Lead_Electrode_001`…, since FreeCAD strips trailing digits before
appending its counter). `tissue_map.yaml` still matches them; `leads.txt` is
what tells you which contact belongs to which lead.

The document is 25 MB, so ten configurations must not mean ten documents, and
every save of one is another full copy in git history. The durable artefact is
therefore the YAML entry plus its STL directory — a dozen numbers and a few
hundred kB. Bodies in a FreeCAD document are a **disposable preview**: the macro
(or `build_lead_config.preview_set()`) drops every lead of the set into the
active document inside a group named `SCS_LeadPreview`, replacing any previous
preview, and “Remove preview” puts the tree back. Explore, export the one you
want, close the document without saving.

Preview `Name`s are prefixed per lead — `SCS_Preview_Lead2_Contact_03` — and that
is not cosmetic. FreeCAD uniquifies a colliding `Name` by **stripping its
trailing digits** and appending a counter, so a second lead's
`SCS_Preview_Contact1` would come back as `SCS_Preview_Contact001` and the two
leads' contacts would alphabetise into each other. That is exactly what once made
one 8-contact lead read as a 12-contact one in the tree.

## What happened to RADO's own 4-contact DRG lead

It stays in `STL_files/` — five files RADO ships, unmodified and tracked — and it
is **removed from the anatomy document** by `drop_rado_lead.py`. The reasoning is
the same one that deleted it from the old study variants, now applied to the
document those were copies of: Simpleware or Gmsh will mesh it, and the solve
will treat four platinum cylinders and a sheath as conductors sitting in the
field a few millimetres off the cord. Hiding a body changes none of that.

Nothing is lost, because the geometry was never *in* the document in any sense
other than “imported from those STLs”. It comes back into any open document in
one click — `build_lead_config.preview_rado_lead()`, or the panel's **Show RADO's
lead** button — as its own disposable overlay group, imported with no transform,
so a generated DRG lead can be looked at right beside it. And it is reproduced as
parameters by `rado_drg_L3_match` in the catalogue, which is what `--compare-rado`
measures against.

**The tradeoff, stated plainly:** `NBF_RADO-SCS.FCStd` stops being a faithful
mirror of what RADO published — someone opening it will not see the lead the
paper's figures show, and has to know about the button. Against that, the
document becomes usable as the mesh/solve input it is meant to be without a step
somebody has to remember, and the step people forget is the one that silently
corrupts a solve rather than failing it.

`drop_rado_lead.py` **saves the document**, so it refuses to run unless
`FreeCAD.GuiUp` — see below. Run it from FreeCAD's Python console, or
`--dry-run` it anywhere.

**Nothing in this pipeline saves a model document headlessly.** A `.FCStd`
written with no GUI running has no `GuiDocument.xml` and no `ShapeAppearance`
blobs at all — every tissue colour gone, and the document still opens looking
like a fresh grey import. `make_tissue_groups.py` documents the measurement;
`build_lead_config.preview_set()` and `preview_rado_lead()` refuse outright to
touch a saved document when `FreeCAD.GuiUp` is false, and `drop_rado_lead.py`,
which is the one script here that does save, refuses for the same reason.

## The macros

Two GUI panels. `lead_designer` is the clinical front end -- pick or dial in a
lead configuration, add up to `max_leads_per_type` leads of each type, check that
each fits *and* that no two of them collide, preview them all in the open model,
and overlay RADO's own lead to compare against.
`tissue_visibility` gives per-tissue show / hide / isolate, which is what makes
the model legible: hiding bone, discs and vessels turns a dense mess into a
clear view of the lead in the epidural space.

**Install -- the macro directory is VERSIONED.** On FreeCAD 26.3 it is
`~/.local/share/FreeCAD/v26-3/Macro/`. A file dropped in
`~/.local/share/FreeCAD/Macro/` is ignored completely and the macro simply never
appears in the list, with no error to tell you why. The path moves with the
FreeCAD version, so ask rather than hardcode:

    MDIR=$(python3 -c 'import FreeCAD, os; print(os.path.join(FreeCAD.getUserAppDataDir(), "Macro"))')
    cp src/freecad/*.FCMacro "$MDIR"

Then in FreeCAD: **Macro → Macros…**, select it, **Execute**. That is all that is
needed to run one.

**Putting one on a toolbar** is a separate, optional step, and it only works
after the macro is in the directory above: **Tools → Customize → Macros** to
give it an icon, then **Tools → Customize → Toolbars** to create a toolbar in
whichever workbench you use and move the macro into it. Toolbar changes take
effect the next time that workbench loads.
