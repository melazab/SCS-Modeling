"""Shared configuration for the FreeCAD-side (gmsh + P1 FEM) SCS pipeline."""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STL = os.path.join(ROOT, "STL_files")
LEAD = os.path.join(ROOT, "fem", "leads", "fem_dorsal_T10")
OUT = os.path.join(ROOT, "fem", "out")

# Tissue bodies of the stripped model, OUTERMOST FIRST in file order but the
# classification order below is INNERMOST FIRST (first match wins).
TISSUE_STL = {
    "epidural": os.path.join(STL, "T8-10 - neuro_EpiduralSpace-1.STL"),
    "dura":     os.path.join(STL, "T8-10 - neuro_Meninges-1.STL"),
    "csf":      os.path.join(STL, "T8-10 - neuro_CSF-1.STL"),
    "white":    os.path.join(STL, "T8-10 - neuro_whitemater-1.STL"),
    "grey":     os.path.join(STL, "T8-10 - neuro_GreyMater-1.STL"),
}
CONTACT_STL = {i: os.path.join(LEAD, f"SCS Lead Electrode {i}.stl") for i in range(1, 9)}
INSULATOR_STL = os.path.join(LEAD, "SCS Lead Insulator.stl")

# Priority for point classification: first hit wins. Metal before insulator
# (contacts are recessed in the insulator body), then inner tissues outward, so
# the 0.06 % of volume where neighbouring shells overlap by a tessellation
# sliver is resolved deterministically.
ORDER = (["contact%d" % i for i in range(1, 9)]
         + ["insulator", "grey", "white", "csf", "dura", "epidural"])

# Electrical conductivity, S/m. EVERY value is taken verbatim from
# src/ansys/tissue_map.yaml (sigma_S_per_m); nothing here is invented.
SIGMA = {
    "grey": 0.276,       # grey_matter
    "white": 0.1432,     # white_matter, ISOTROPIC as in tissue_map (see README)
    "csf": 1.7,          # csf
    "dura": 0.037,       # meninges_dura
    "epidural": 0.04,    # epidural_space
    "insulator": 2.0e-5,  # lead_insulation
    "background": 0.04,  # tissue_map "vertebra" -- only used in the variant run
}
# tissue_map gives metal electrode 4.0e6 S/m. A 1e8 conductivity ratio against
# epidural fat wrecks the conditioning of the linear system for no physical
# gain: at 1e4 S/m a contact is already equipotential to ~1e-5 of the driving
# voltage. SIGMA_METAL_TRUE is kept so the clamp can be audited.
SIGMA_METAL_TRUE = 4.0e6
SIGMA_METAL = 1.0e4

# Bipolar drive, as specified: source contact, sink contact, current in amperes.
SOURCE_CONTACT = 3
SINK_CONTACT = 5
CURRENT_A = 1.0
