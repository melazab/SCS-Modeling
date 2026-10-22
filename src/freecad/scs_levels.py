"""Where each vertebra sits along the cord: one definition for every SCS panel
and script.

A level's position is its vertebra's VOLUME centroid projected square onto the
cord centreline (fem/out/geometry/rado_up4/centreline.json), as arc length in
mm from the centreline's caudal end. On the extended model (2026-10-22):

    T6 212.0  T7 181.0  T8 150.5  T9 119.5  T10 89.0  T11 58.0  T12 26.0

evenly spaced at the build's 30.7 mm level step. The vertebra labels
(add_vertebra_labels.py) point at the right transverse process tip instead,
2-6.5 mm cranial of these and not consistently so (T10 -> T9 reads 35 mm), so
they are annotations, not positions. Before this module the panels projected
the labels and the tract-body build projected their z only, and the two
disagreed by up to a millimetre; the Cross-section viewer's thin slice then
missed the Tract visualizer's per-vertebra slice at T8 (Mohamed, 2026-10-22).

Pure numpy; FreeCAD is only needed for document_levels().
"""
import re

import numpy as np

# Extended-set vertebra STLs by level (src/freecad/body_aliases.yaml names
# RADO's three; extend_rado.py's copies of V1-4 are __up1.. upward).
VERTEBRA_STEMS = {"T12": "V1-2", "T11": "V1-3", "T10": "V1-4", "T9": "V1-4__up1",
                  "T8": "V1-4__up2", "T7": "V1-4__up3", "T6": "V1-4__up4"}
LABEL_RE = re.compile(r"^Vertebra (T\d+)$")


def volume_centroid(verts, tris):
    """Centroid of the solid a closed triangle mesh bounds."""
    v = np.asarray(verts, float)
    t = np.asarray(tris, int)
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    vol = np.einsum("ij,ij->i", a, np.cross(b, c))/6.0
    return ((a + b + c)/4.0*vol[:, None]).sum(0)/vol.sum()


def arc_length(point, s_mm, points_mm):
    """Arc length of the centreline sample nearest to point in the sagittal
    (y, z) plane. The centreline lies in x = x0, so this is the foot of the
    perpendicular from the point to the cord, to the sampling (0.5 mm)."""
    P = np.asarray(points_mm, float)
    k = int(np.argmin((P[:, 1] - point[1])**2 + (P[:, 2] - point[2])**2))
    return float(np.asarray(s_mm, float)[k])


def stl_levels(geometry_dir, s_mm, points_mm, read_stl):
    """Levels from the extended set's vertebra STLs (no FreeCAD)."""
    import glob
    import os
    out = {}
    for level, stem in VERTEBRA_STEMS.items():
        hits = [p for p in glob.glob(os.path.join(geometry_dir, "*.STL"))
                if os.path.basename(p).split(" - ", 1)[-1].rsplit("__conformal", 1)[0] == stem]
        if len(hits) == 1:
            out[level] = round(arc_length(volume_centroid(*read_stl(hits[0])), s_mm, points_mm), 2)
    return out


def document_levels(doc, s_mm, points_mm):
    """Levels from the open document's "Vertebra Tn" bodies."""
    out = {}
    for o in doc.Objects:
        m = LABEL_RE.match(o.Label)
        if not m:
            continue
        if hasattr(o, "Mesh"):
            cen = volume_centroid(*o.Mesh.Topology)
        elif hasattr(o, "Shape") and not o.Shape.isNull():
            cen = np.array(o.Shape.CenterOfMass)
        else:
            continue
        out[m.group(1)] = round(arc_length(cen, s_mm, points_mm), 2)
    return out
