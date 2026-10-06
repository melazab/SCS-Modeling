# Reverse-engineering pilot — Vertebra 01

Source: `T8-10 - V1-2.STL`. Work is scoped to this one vertebra at Mohamed's
request; the full inventory and other pilots remain deferred. Assemblies of
parts are permitted. Anatomical level remains unassigned because source and
repository level labels conflict. Source files are unchanged.

**Current candidate does not meet acceptance.** The editable vertebral body
has an 18-pole periodic B-spline sketch → Pad → rim Fillet → two-channel Pocket.
The posterior is a separate static faceted solid, not reconstructed parametric
anatomy. Both components are valid closed solids, but overlap; native CAD
fusion failed. No conformal shared interface is claimed.

- [Native candidate assembly](pilot/vertebra_01_candidate/Vertebra_01_Assembly_Candidate.FCStd)
- [Editable body alone](pilot/vertebra_01_candidate/Vertebra_01_Body_Candidate.FCStd)
- [Detailed measured report and limitations](pilot/vertebra_01_candidate/REPORT.md)
- [Candidate preview](pilot/vertebra_01_candidate/candidate.png)

| Measurement | Result |
|---|---:|
| Diagnostic source → candidate sampled RMS | 0.03066854 mm |
| Diagnostic candidate → source sampled RMS | 0.05416842 mm |
| Equal-direction pooled sampled RMS | 0.04401578 mm |
| Bidirectional sampled maximum | 1.82762902 mm |
| Diagnostic union volume error | -0.123533% |
| Candidate investigation elapsed time | 11.18 minutes |

These distance and volume results are for a **diagnostic mesh union**, not a
successful native CAD fusion. Distances use 100,000 area-uniform samples per
direction; the sampled maximum is not a certified Hausdorff distance. The
large error at the posterior attachment/cut transition rules out accepting
this candidate. The investigation time excludes initial context reading and
the earlier teaching session; it is not a completed-body production estimate.

Editable parameters include pad height, rim radius, channel radii/centers,
spline poles, and Placement. Height +1 mm, rim radius +0.1 mm, and one channel
radius +0.05 mm each passed validity checks and were restored. Anatomical
width/depth controls and posterior adaptation are not implemented.

Next engineering step: fit the body–pedicle transition from multiple source
sections and construct an explicit shared CAD boundary. Do not extend the
current approximation to other vertebrae.

`pilot/vertebra_01/` retains preparation code, measurements, and derived source
references used by the candidate. Educational notes, demo macro, screenshot,
demo results, and lesson FCStd were removed at Mohamed's request. The user's
open FreeCAD session was not modified during cleanup. Engineering scripts,
candidate models, validation data, and reports are retained; nothing committed
or pushed.
