# `src/ansys/` — solver side

The electric-conduction finite-element work: MAPDL decks, the SLURM job trees
they ran under, the tissue/material map, and post-processing of results.

## Install

    pip install -r requirements.txt      # only needed for plot_voltage_slice.py

Ansys itself is licensed software on the CWRU Pioneer cluster:

    module load ansys/25.1               # MAPDL is then `ansys251`

## The tissue map

`tissue_map.yaml` is the single source of truth for what every body is made of:
STL filename patterns → tissue class → Ansys Engineering Data material name →
resistivity/conductivity → display colour. Both halves of `src/` read it.

    python3 check_tissue_map.py          # asserts all 245 bodies map to exactly one tissue

**Never change a conductivity here to make a solve behave.** The values are
matched to Khadka et al. on purpose. The deck is `/units,uMKS`, so resistivity in
it is in Tohm·µm — `MP,RSVX,1,2.5e-13` is 2.5e-7 Ω·m, which *is* the Engineering
Data value for platinum-ish contacts, not a typo. That misreading cost a wasted
overnight run.

## Running a solve

Job directories are self-contained; each holds its own `run_*.sbatch` with a
header explaining what that run was testing and what it found.

    sbatch run_pinball_smp.sbatch

Cluster facts worth knowing before you submit:

- The `tlv` account is capped at **24 CPUs**, and SLURM bills CPUs in proportion
  to memory (~6 GB/core on `batch`, ~29 GB/core on `smp`). One 24-CPU job
  therefore blocks every other job you have, including your own interactive
  `srun`.
- `batch` tops out around 140 GB for that cap; `smp` nodes carry ~1.1 TB. The
  full 20.8M-node model needs `smp` for the direct sparse solver.
- Plain `-np N` **fails** on this install. Use `-smp` for small jobs and
  `-dis -mpi intelmpi -np N` for large ones.
- Use a login shell so `module` exists: `ssh case-hpc 'bash -ls' <<'EOF' … EOF`.

## Deck gotchas

The production deck is 4.5 GB — never `cat` it.

- It has **CRLF** line endings, so `sed` patterns anchored with `$` silently
  never match while unanchored ones do. Strip `\r` first and *assert* that every
  patch landed; the sbatch scripts here all do, and abort if one didn't.
- `/units,uMKS` means ±1 mA is written as ±1e9 (picoamps).
- `EDELE` is invalid inside `/SOLU` and is silently ignored.
- `*GET,…,MXV` and `*GET,…,cp,0,maxi` are invalid labels that abort `/POST1` in
  batch mode. Use `NSORT` then `*GET,…,SORT,0,MAX`.

## Headless Mechanical

Works, but its scripting engine is **IronPython 2.7** — Python 3 syntax fails
silently with zero output and exit code 0:

    /usr/local/ansys_inc/v251/aisol/.workbench -DSApplet -AppModeMech -b -script /abs/path.py

Messages are at `ExtAPI.Application.Messages`; nothing useful reaches stdout.

## Anisotropic white matter

Engineering Data holds white matter as **isotropic** 0.1432 S/m (the paper's
transverse value). The paper's anisotropy — 0.1432 transverse / 0.6 longitudinal
— cannot be entered through the Mechanical Electric GUI. It needs an APDL command
snippet (`MP,RSVX/RSVY/RSVZ`) or PyMAPDL; see `mapdl_electric_solve_template.dat`.
The cord axis is global Z in this STL set. Currently left isotropic by choice.
