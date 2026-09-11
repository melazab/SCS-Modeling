# Pre-flight — validate a patched deck before spending an hour on it

Three full-model runs were lost to bugs that were each detectable in seconds:

| job | bug | cost |
|---|---|---|
| 3797286 | `--mem` below actual need | 2h35m, OOM-killed |
| 3802282 | `finish` dropped to BEGIN, so `solve` was ignored | 5:48, no solve at all |
| 3804186 | `CPNGEN` issued in `/SOLU`, sets came out 1 node each | 35:29, wrong model solved |

They are the same class of mistake. Every patch was verified by grepping the
text — "did the sed land?" — and none was verified to be *valid APDL in the
right processor*. A `grep -q` cannot tell you that `CPNGEN` is a `/PREP7`
command.

**Rule: nothing goes to the full model until it has run here first.**

## `testA` is the harness, and it already works

`../testA_mapdl_cylinder/testA.dat` is a small coaxial cord model -- cord in
CSF in fat, a 1.25 mm lead with two 3 mm contacts -- and it **solves**:
40,274 nodes in the exported slab, voltage −1.0366 to +1.0380 V for a ±1 mA
bipolar pair. It exercises exactly the patterns the real model needs:
SOLID232, anisotropic white matter, equipotential contacts, current injection,
a reference potential.

It runs in minutes on a couple of cores. Any command sequence intended for the
4.5 GB deck should be rehearsed here, where being wrong costs nothing.

## The idiom that works, and the one that does not

testA has made an equipotential electrode correctly all along:

    CMSEL,S,CONTACT1
    NSLE,S
    CP,NEXT,VOLT,ALL          ! one set, every selected node
    *GET,NA,NODE,0,NUM,MIN
    F,NA,AMPS,1.0e-3

One `CP` over a *selected node set*. What job 3804186 did instead was emit one
`CPNGEN` per node, which produced 878 coupled sets of one node each:

    COUPLED SET=     1  DIRECTION= VOLT  TOTAL NODES=     1
    MAXIMUM COUPLED SET NUMBER=     0

Worse than useless, because the conversion had already deleted the 878 `CE`
lines, so that run had no equipotential constraint on the electrodes at all.
It only looked plausible because the meshed platinum enforces equipotential by
itself -- the very redundancy the change existed to remove.

So the conversion should select the electrode's nodes and issue a single `CP`:

    nsel,none
    nsel,a,node,,786
    nsel,a,node,,787
    ...
    cp,next,volt,all

And the job must assert `TOTAL NODES` is the expected count, not merely that
some `CP` appeared.

## What testA also tells us about the real problem

Read its own header comment:

    moderated electrode conductivity (1e4 S/m instead of 4e6) and NO modeled
    insulator body -> avoids the 1e11 material contrast that produced the
    negative-pivot failure in the 20.8M-node runs

testA works *because* the conductivity contrast was taken out. It is a
standing demonstration that the physics recipe is sound and the contrast is
what breaks it -- which is the same conclusion the ECC sweep and the stripped
run reached independently, from the other direction.

Note also what testA does NOT have: any contact pairs. Its nested cylinders are
`VOVLAP`'d and `NUMMRG`'d into a single conformal mesh with shared nodes at
every interface. No CONTA174, no pinball, no ECC, no constraint equations. That
is how the published workflow meshes these models (Simpleware for Khadka et
al.), and it is the structural difference between testA, which solves, and the
production deck, which does not.
