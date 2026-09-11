# ECC sweep (job 3800027) — ECC is not the lever

All three array tasks **COMPLETED, exit 0:0**, ~40–50 min each, ~420 GB.
The patch definitely applied — each run logged
`ECC blocks set to 4.e+0NN: 1116`, `pinball pairs patched: 2232`, and
`materials confirmed at original Khadka values`.

| `_maxCond` | max pivot | min pivot | result |
|---|---|---|---|
| 4e+12 | 5.555441222E+15 | −4.617862376E−06 | terminated |
| 4e+10 | 5.555441222E+15 | −4.617862376E−06 | terminated |
| 4e+08 | 5.555441222E+15 | −4.617862376E−06 | terminated |

**Bit-identical across four orders of magnitude.** Not "similar" — the same
digits. Whatever sets the pivots, it is not the electric contact conductance.
The hypothesis that motivated this sweep is dead, and no amount of further
ECC tuning will change anything.

One thing did change versus the 4e+16 run (job 3797962, max pivot 7.13e16):
the maximum pivot dropped to 5.55e15 and then **stopped moving**. So ECC does
influence the matrix until it falls below whatever now dominates at 5.55e15,
and every value tested here is already below that. The plateau was reached
before the first sample.

## What the failure actually is now

    A large negative pivot value ( -4.617862376E-06 ) has been encountered
    in the global assembled matrix at the VOLT degree of freedom of node
    12254130.

Not the "insufficiently constrained model" of the original runs — that was
the floating-body signature and the pinball fix retired it. This is a
**negative** pivot on a voltage DOF, with a pivot ratio around
5.55e15 / 4.6e-06 ≈ 1.2e21.

## The 878 constraint equations are electrodes, and must not be deleted

    ce,next,0.,3149,volt,-1.,3150,volt,1.
    ce,next,0.,3149,volt,-1.,3151,volt,1.

878 constraint equations, **every one of them on `volt`**, each tying a
contact-surface node to a master node. This is how Workbench makes an
electrode equipotential: one conductor, one potential.

That settles an old loose end. `CEDELE,ALL` once "unblocked" a truncated
solve, and it was recorded as a fix. It is not one — it works by deleting the
electrode definitions. Every contact stops being equipotential, so the model
still solves but it is no longer the model. **Do not use `CEDELE,ALL`.**

It also explains the sign. Equipotential constraints enter as Lagrange
multipliers, which make the assembled system a saddle-point problem rather
than positive-definite, and negative pivots are a normal feature of those.
MAPDL tolerates them until the ratio gets extreme. Ours is 1e21.

## Contact settings, for the record

    keyo,cid,1,6     DOF set (electric)
    keyo,cid,2,0     augmented Lagrange
    keyo,cid,12,5    bonded, always
    rmod,cid,6,-500  pinball, absolute 500 um   (our change)
    rmod,cid,19,_maxCond/_ASMDIAG               (ECC — no longer believed to matter)

Note ECC is set as `_maxCond/_ASMDIAG`, a ratio, so the effective value
depends on an assembled diagonal we do not control. That is a plausible
reason the absolute `_maxCond` we set stopped mattering.

## Where the conditioning actually comes from

The prime remaining suspect is the conductivity range the model spans, which
is enormous and entirely physical:

    metal electrode   4e6    S/m
    lead insulation   2e-5   S/m      -> a ratio of 2e11
    CSF               1.7    S/m
    white matter      0.1432 S/m

**These values are correct and are not to be changed.** They are matched to
Khadka et al.; that is settled and was settled for good reason.

The way out is not to alter the physics but to stop meshing the part that
causes the trouble. A 4e6 S/m platinum cylinder is, electrically, an
equipotential surface — which is exactly what the 878 constraint equations
already say. Modelling it as a *boundary condition* on the contact surface,
rather than as a meshed volume of extremely conductive material, is standard
practice in bioelectric FEM and removes eleven orders of contrast from the
matrix without changing what is being simulated.

That is a larger change than a `sed` on the deck: it means rebuilding the
model with the lead as surfaces rather than solids. It also happens to suit
where the project is going, since the leads are now parametric and generated
per configuration rather than baked into the geometry.

## Recommended next step

Do not submit another full-model solve yet. Test the boundary-condition idea
on the truncated sub-model first, where a run costs minutes instead of an
hour and 420 GB: replace the meshed electrode volumes with surface boundary
conditions, confirm the negative pivot goes away and the pivot ratio comes
back to something sane, and only then take it to the full model.
