"""Recover Lead Designer parameters from the leads already in a document.

WHY THIS EXISTS
---------------
Lead Designer does not persist its parameters. A built lead becomes plain
`Part::Feature` shapes named `SCS_Preview_LeadN_Contact_MM` and
`SCS_Preview_LeadN_Insulator`, carrying only Placement, Shape, Label and
Visibility -- verified against a live document; there is no custom property
anywhere holding contacts/length/gap/diameter. So opening a RADO model that
already carries leads showed the panel's DEFAULTS, which described nothing in
front of you: the poster dorsal model has TWO leads at +/-2.00 mm from the
midline, 1.60 mm across, centred at z = 80.44, and the panel offered one lead
on the midline, 1.30 mm across, at z = 110.

This module reads the geometry back instead.

THE FRAME, AND WHY THE NAIVE MEASUREMENT IS WRONG
-------------------------------------------------
The designer lays contacts out along **z**: centre-to-centre z pitch is exactly
`contact_length + gap` (measured 4.0000 mm on the poster dorsal lead). But the
lead is swept along a CURVED centreline, so anything measured along the lead's
own local axis is longer than its z counterpart by the same factor -- 1.038 for
that lead. Measuring a contact's length along the local axis and reading it as
the designer's `contact_length` therefore overstates it by ~4% (3.116 vs 3.0),
and the same for the gap.

The fix is to measure the RATIO rather than the length. `axis_length /
axis_pitch` cancels the curvature factor, and multiplying by the exactly-known
z pitch returns the designer's own numbers: **3.002 mm and 0.998 mm against a
true 3.0 and 1.0.** That is the whole trick, and it is why this does not simply
read bounding boxes.

WHAT IS AND IS NOT RECOVERABLE
------------------------------
  contacts, contact_length, gap, diameter, x_offset, z_center   recovered
  tail                           recovered from the insulator's overhang
  type (dorsal / ventral / drg)  read from the object Label
  contact_thickness, level, note NOT encoded anywhere -- defaults are kept

`type` comes from the Label because nothing else in the file records it. Labels
can be stale: the poster VENTRAL document still labels its contacts
"lead 1 (dorsal)", inherited from the Designer. `recover()` therefore reports
the label it used in each result's `_label_type` so a caller can flag it, and a
geometric cross-check lives in `infer_side_from_cord()`.
"""
import math
import re

_CONTACT_RE = re.compile(r"^SCS_Preview_Lead(\d+)_Contact_(\d+)$")
_INSULATOR_RE = re.compile(r"^SCS_Preview_Lead(\d+)_Insulator$")
_LABEL_TYPE_RE = re.compile(r"lead\s+\d+\s*\((\w+)\)", re.I)

# build_lead_config.MIDLINE_X. Duplicated rather than imported so this module
# stays importable without the builder; recover() takes an override.
MIDLINE_X = 56.60

TESSELLATION_MM = 0.05


def _by_lead(doc):
    contacts, insulators = {}, {}
    for obj in doc.Objects:
        m = _CONTACT_RE.match(obj.Name)
        if m:
            contacts.setdefault(int(m.group(1)), {})[int(m.group(2))] = obj
            continue
        m = _INSULATOR_RE.match(obj.Name)
        if m:
            insulators[int(m.group(1))] = obj
    return contacts, insulators


def _points(obj):
    """Surface points of a shape, via tessellation.

    A swept contact has very few topological Vertexes (a cylinder has two), so
    measuring from `Shape.Vertexes` would sample almost nothing. Tessellating
    gives the surface.
    """
    return obj.Shape.tessellate(TESSELLATION_MM)[0]


def _axis_metrics(centres, obj, k):
    """(length along the local axis, diameter) for contact k."""
    a = centres[max(k - 1, 0)]
    b = centres[min(k + 1, len(centres) - 1)]
    axis = b.sub(a)
    if axis.Length < 1e-9:                       # single contact: fall back to z
        import FreeCAD
        axis = FreeCAD.Vector(0.0, 0.0, 1.0)
    axis.normalize()
    along, radial = [], []
    for p in _points(obj):
        d = p.sub(centres[k])
        t = d.dot(axis)
        along.append(t)
        radial.append(math.sqrt(max(d.Length ** 2 - t * t, 0.0)))
    return max(along) - min(along), 2.0 * max(radial)


def recover(doc, midline_x=MIDLINE_X):
    """Parameters for every Lead Designer lead in `doc`, lowest lead number first.

    Returns a list of dicts shaped like build_lead_config's parameter sets, each
    with the extra keys `_lead_number`, `_label_type` and `_recovered` (the list
    of parameter names that came from geometry rather than from a default).
    Returns [] when the document holds no Lead Designer leads.
    """
    if doc is None:
        return []
    contacts, insulators = _by_lead(doc)
    out = []
    for n in sorted(contacts):
        idx = sorted(contacts[n])
        if not idx or idx != list(range(1, len(idx) + 1)):
            continue                              # not a complete lead; skip it
        objs = [contacts[n][i] for i in idx]
        centres = [o.Shape.CenterOfMass for o in objs]

        lengths, diameters = [], []
        for k, o in enumerate(objs):
            L, d = _axis_metrics(centres, o, k)
            lengths.append(L)
            diameters.append(d)
        axis_len = sum(lengths) / len(lengths)
        diameter = sum(diameters) / len(diameters)

        if len(centres) > 1:
            axis_pitch = sum(centres[k + 1].sub(centres[k]).Length
                             for k in range(len(centres) - 1)) / (len(centres) - 1)
            z_pitch = sum(centres[k + 1].z - centres[k].z
                          for k in range(len(centres) - 1)) / (len(centres) - 1)
            # The ratio cancels the curvature factor; the z pitch is exact.
            frac = axis_len / axis_pitch if axis_pitch > 1e-9 else 0.75
            contact_length = frac * abs(z_pitch)
            gap = abs(z_pitch) - contact_length
        else:
            contact_length, gap = axis_len, 0.0

        xs = [c.x for c in centres]
        zs = [c.z for c in centres]
        params = {
            "contacts": len(idx),
            "contact_length": round(contact_length, 3),
            "gap": round(gap, 3),
            "diameter": round(diameter, 3),
            "x_offset": round(sum(xs) / len(xs) - midline_x, 3),
            "z_center": round(sum(zs) / len(zs), 3),
        }
        recovered = sorted(params)

        ins = insulators.get(n)
        if ins is not None:
            # tail is plain insulator beyond the END CONTACTS, per end
            # (build_lead_config: total = array_length + 2 * tail).
            ins_lo, ins_hi = ins.Shape.BoundBox.ZMin, ins.Shape.BoundBox.ZMax
            arr_lo = min(o.Shape.BoundBox.ZMin for o in objs)
            arr_hi = max(o.Shape.BoundBox.ZMax for o in objs)
            tail = ((arr_lo - ins_lo) + (ins_hi - arr_hi)) / 2.0
            if tail > 0:
                params["tail"] = round(tail, 3)
                recovered.append("tail")

        m = _LABEL_TYPE_RE.search(getattr(objs[0], "Label", "") or "")
        label_type = (m.group(1).lower() if m else None)
        params["type"] = label_type if label_type in ("dorsal", "ventral", "drg") else "dorsal"
        params["_lead_number"] = n
        params["_label_type"] = label_type
        params["_recovered"] = recovered
        out.append(params)
    return out


def infer_side_from_cord(doc, params):
    """Geometric cross-check on `type`, for when the Label is stale.

    Returns "dorsal", "ventral" or None. The cord's own dorsoventral axis is
    taken from the white-matter mesh at the lead's z: contacts at higher y than
    the cord are dorsal to it. Returns None when no white-matter object is
    found, which is the common case for a document holding only leads.
    """
    white = None
    for obj in doc.Objects:
        if "whitemater" in obj.Name.lower() or "whitematter" in obj.Name.lower():
            white = obj
            break
    if white is None or not hasattr(white, "Mesh"):
        return None
    z = float(params["z_center"])
    ys = [p.y for p in white.Mesh.Points if abs(p.z - z) < 1.0]
    if not ys:
        return None
    cord_y = sum(ys) / len(ys)
    lead_y = None
    contacts, _ = _by_lead(doc)
    objs = contacts.get(params.get("_lead_number"))
    if objs:
        cs = [objs[i].Shape.CenterOfMass.y for i in sorted(objs)]
        lead_y = sum(cs) / len(cs)
    if lead_y is None:
        return None
    return "dorsal" if lead_y > cord_y else "ventral"
