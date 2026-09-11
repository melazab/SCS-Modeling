# Stripped run (job 3802339) — it is conditioning, not size

Vertebrae, vasculature and intervertebral discs removed by named component.
The strip ran correctly (`Assembly VERTEBRA / VASCULATURE /
INTERVERTEBRAL_DISK is deleted`) and the solve ran to factorization.

|                     | full (3797962) | stripped (3802339) |
|---------------------|----------------|--------------------|
| elements            | 12,889,916     | **11,397,737**     |
| nodes               | —              | 8,621,408          |
| peak memory         | 419 GB         | **263 GB**         |
| elapsed             | 37 min         | **26 min**         |
| max pivot           | 7.13e16        | **7.73e16**        |
| min pivot           | −1.56          | −7.96e−06          |
| outcome             | terminated     | terminated         |

**The pivot ratio did not improve.** A third less memory, 1.5 million fewer
elements, and the maximum pivot went slightly *up*. That settles the question
the run was submitted to answer: the failure is conditioning, not problem
size, and no amount of stripping anatomy will fix it.

The strip is still worth keeping for what it was really for — 26 min and
263 GB instead of 37 min and 419 GB makes the next experiments affordable.

Note the element count fell only 11.6%, not the ~37% a material-based
estimate suggested. That estimate was wrong for the reason the run was
careful to avoid: material 6 covers vertebrae **and** epidural space
together, so counting mat-6 elements as "vertebrae" overstated it.

## The smoking gun

    Distributed sparse solver maximum pivot = 7.733301104E+16 at node 3149

Node **3149** is an electrode constraint-equation master — it appears as the
master in 505 of the deck's 878 `ce,...,volt` lines:

    ce,next,0.,3149,volt,-1.,3150,volt,1.

So the largest pivot in an 11.4-million-element model sits exactly on the
node that holds an electrode at a single potential. Combined with the ECC
sweep (three values four orders apart, bit-identical pivots), the electrode
representation is where the conditioning is, not the tissue.

The failure itself is still a negative pivot on a VOLT degree of freedom,
now at node 12281509, magnitude −7.96e−06.

## Why this supports the boundary-condition plan

The electrode is currently represented **twice**: as a meshed volume at
4e6 S/m, and as constraint equations forcing that volume to one potential.
The second representation is the physically meaningful one. The first
contributes only a 2e11 conductivity ratio against the insulation at
2e-5 S/m.

Two changes follow, and this run supports both:

1. **Stop meshing the metal.** Impose each contact as a surface condition
   instead. Removes the 2e11 contrast that drives the maximum pivot.
2. **Use `CP`, not `CE`.** All 878 constraint equations are the trivial form
   "node B has the same VOLT as node A" — which is exactly `CP,,VOLT`. The
   deck contains **878 `ce,` lines and zero `cp,` lines**. `CE` is a general
   linear constraint and can bring Lagrange multipliers, making the system
   indefinite, which is where negative pivots come from. `CP` is implemented
   by elimination: coupled DOFs collapse into one unknown, so the matrix gets
   *smaller* rather than harder.

Passive contacts stay. A floating conductor is real and shunts current — it
is Khadka's own boundary condition ("floating on the rest"). It becomes a
`CP` set with no current injected, which costs less than the meshed platinum
plus constraint equations it replaces.

Unverified and worth checking on the truncated sub-model first: that MAPDL
eliminates `CP` rather than falling back to a multiplier for the VOLT DOF in
this configuration.
