# What SCS actually stimulates — a primer

Written to fill in the neuroscience big picture before the NEURON work in
`neuron_plan.md`. Everything here is standard, textbook-level material; the
references at the bottom are the ones to cite in a poster or paper.

## The one-line answer

**Axons. Almost never cell bodies.**

Extracellular stimulation activates *myelinated axons* at far lower amplitudes
than it activates somata or dendrites. So when a clinician turns on an SCS
system, what responds first is a bundle of large myelinated fibres passing by the
electrode. The cell bodies those fibres belong to are usually centimetres away
and are not directly involved at all.

### Why axons, biophysically

A fibre is driven not by the extracellular potential *V<sub>e</sub>* but by its
**second spatial derivative along the fibre** — Rattay's activating function:

        f  ∝  d · ∂²V_e/∂x²

Three consequences that matter for everything below:

1. **Bigger fibres are recruited first.** Threshold falls roughly with the
   inverse of diameter. This inverts the physiological recruitment order — normal
   sensory input recruits small fibres first, SCS recruits the largest first.
2. **Orientation matters as much as distance.** A fibre bending through a field
   gradient sees a large ∂²V<sub>e</sub>/∂x² and is easy to recruit; a straight
   fibre in a locally uniform field is hard. This single fact explains most of
   what follows.
3. **Cathodic stimulation depolarises the node under the electrode** and
   hyperpolarises the flanking nodes, which is why cathodes are the working
   contacts.

## The dorsal target

### What is in the dorsal columns

The dorsal columns are the ascending white-matter tracts carrying **touch,
vibration and proprioception**. The fibres in them are the **central processes of
primary sensory neurons**, and this is the part worth being precise about:

- The cell body sits in the **dorsal root ganglion**, *outside* the spinal cord.
- It is **pseudounipolar**: one process goes to the periphery, one enters the
  cord through the dorsal root.
- Once inside, that central process **bifurcates** — one branch ascends in the
  dorsal column toward the brainstem, one descends, and **collaterals** dive into
  the dorsal horn to synapse locally.

So conventional SCS stimulates a fibre whose soma is in the DRG and whose synapse
is in the dorsal horn, at a point on the axon in between. The soma is not in the
loop.

At **T8–T10** — the levels in this model, and the standard clinical placement for
low-back and leg pain — the dorsal column here is essentially the **gracile
fasciculus**, the medial tract carrying lower-limb and lower-trunk afferents.
(The cuneate fasciculus, for the arms, only appears above ~T6.) That is precisely
why this lead position treats leg and low-back pain: it is where those fibres run.

### How stimulating touch fibres relieves pain

It does not, directly. There are **two hops**:

1. SCS depolarises large **Aβ** fibres in the dorsal column.
2. Their **collaterals** excite inhibitory interneurons in the dorsal horn, which
   suppress transmission from small nociceptive **Aδ and C** fibres onto the
   projection neurons that carry pain signals to the brain.

This is **gate control** (Melzack & Wall, 1965) — the reason SCS exists. The
stimulation target and the therapeutic site are different structures, one synapse
apart. The paraesthesia the patient feels is the direct effect; the analgesia is
the indirect one.

Worth a caveat in any write-up: gate control describes **conventional,
paraesthesia-based SCS at ~40–60 Hz**. Newer paradigms — 10 kHz, burst — relieve
pain *without* paraesthesia, and their mechanisms are still debated. Do not
present gate control as an explanation for those.

## The real modelling problem: competing populations

If dorsal column recruitment were all that happened, you would simply turn the
amplitude up. The constraint is that **something else is easier to recruit.**

| population | where | recruited | effect |
|---|---|---|---|
| **Dorsal root (DR) fibres** | in the rootlets, before entering the cord | **first — lowest threshold** | sharp segmental/radicular paraesthesia, dermatomal banding, discomfort |
| **Dorsal column (DC) fibres** | in the cord, running longitudinally | second | the *wanted* effect — paraesthesia over the painful area |

Dorsal root fibres win because they are **closer to the electrode**, they **curve**
through the field (large ∂²V<sub>e</sub>/∂x²), and they sit in **highly conductive
CSF**. Dorsal column fibres are further away, run parallel to the electrode array
in a comparatively uniform field, and are shielded behind the CSF layer.

The ratio between these two thresholds is the **therapeutic window** (Holsheimer's
"usage range"):

        therapeutic window  =  DR threshold / DC threshold

Widening it is the entire objective of computational SCS lead design, and it is
the number this project should be reporting. **This is what a recruitment curve is
for**: not "how many fibres fire", but "how much DC recruitment can I buy before
DR recruitment makes it unpleasant".

## What sits between the electrode and the target

Measured by ray-casting this model at midline, z = 110 mm, from the lead inward:

        epidural fat   2.32 mm      lead sits in here
        dura           0.49 mm
        CSF            2.76 mm      <-- the one that matters
        white matter   2.35 mm      dorsal columns: the target
        grey matter    0.44 mm

**The dorsal CSF layer is the single dominant anatomical determinant of SCS
thresholds** (Holsheimer, 2002). CSF is by far the most conductive tissue in the
model — 1.7 S/m against 0.1432 S/m for white matter, more than 10× — so it acts
as a shunt, spreading current longitudinally instead of letting it pass into the
cord. Thicker dCSF means higher thresholds and a narrower therapeutic window. It
varies with spinal level, posture, and between patients, which is a large part of
why SCS outcomes vary. At 2.76 mm this model sits in a realistic human range.

## The ventral side — why it is a different question

Ventrally the geometry is tighter (**1.71 mm** of epidural fat against 2.32 mm
dorsally), and what is nearby is different in kind:

- **Ventral roots** — motor efferents, the axons of **α-motor neurons** whose cell
  bodies *are* in the cord, in the ventral horn. These are among the largest
  fibres in the body (Aα, 13–20 µm), so by the diameter rule they have **very low
  thresholds**. Recruiting them causes muscle contraction, which is the obvious
  risk of ventral stimulation.
- **Anterolateral / spinothalamic tract** — the ascending pathway that actually
  carries nociception. Interesting as a target precisely because it is the
  pain-signalling pathway itself rather than a gate-control proxy.
- **Descending motor tracts** — corticospinal and others in the ventral funiculus.

Ventral epidural stimulation is **not standard clinical practice for pain**, which
is exactly what makes it a research question. The related literature to know is
epidural stimulation for **motor restoration after spinal cord injury**, where
ventral/motor recruitment is the goal rather than the side effect.

For this project the honest framing is: dorsal placement has a well-established
target and a well-defined therapeutic window; ventral placement is being
characterised. The model's job is to say what ventral placement recruits, at what
amplitudes, relative to dorsal.

## Fibre types

| class | diameter | conduction | role | relevance to SCS |
|---|---|---|---|---|
| **Aα** | 13–20 µm | 80–120 m/s | motor, proprioception | lowest threshold; the ventral-root risk |
| **Aβ** | 6–12 µm | 35–75 m/s | touch, vibration | **the target** |
| **Aδ** | 1–5 µm | 5–35 m/s | fast/sharp pain, temperature | thresholds too high to recruit directly |
| **C** | 0.2–1.5 µm | 0.5–2 m/s | slow pain | unmyelinated; effectively never recruited |

SCS works on the **top of this table**, and relieves pain carried by the
**bottom** of it — indirectly, through the dorsal horn.

## What this means for the model

The populations worth simulating, and why:

1. **Dorsal column Aβ fibres** — a range of diameters and depths in dorsal white
   matter. The wanted effect.
2. **Dorsal root Aβ fibres** — following the root trajectories. The constraint
   that sets the upper limit on amplitude.
3. **Ventral root Aα fibres** — the motor side, for the ventral arm of the study.

The headline output is not a single threshold but the **ratio** between (1) and
(2) — the therapeutic window — computed for dorsal and for ventral placement at
each rostro-caudal level.

And the modelling choice this justifies: the **MRG** axon model. These are
myelinated mammalian fibres, so the model needs explicit nodes of Ranvier, myelin
and a double-cable structure. Classical Hodgkin–Huxley is an *unmyelinated squid*
axon at 6.3 °C and is wrong on every one of those counts.

## References worth citing

- Melzack R, Wall PD (1965). *Pain mechanisms: a new theory.* Science — gate control.
- Rattay F (1986/1989). Activating function; why the second derivative of V<sub>e</sub> drives excitation.
- Holsheimer J (2002). *Which neuronal elements are activated directly by spinal cord stimulation?* Neuromodulation — the dCSF result and the DR/DC window.
- Struijk JJ, Holsheimer J, Boom HBK (1993). Recruitment of dorsal column fibres in SCS.
- McIntyre CC, Richardson AG, Grill WM (2002). *J Neurophysiol* — the MRG axon model. ModelDB #3810.
- Lempka SF, Patil PG (2018). Innovations in SCS — a good modern review.
- Khadka N et al. — the RADO-SCS model this geometry comes from.

Verify page numbers and exact titles before they go into a bibliography.
