"""Label each vertebra (T6-T12, or whatever levels the document holds) in the 3D view.

Run it from the FreeCAD GUI with the anatomy document active: Macro > Macros...,
or in the Python console

    exec(open('/home/mohamed/Projects/SCS-Modeling/src/freecad/add_vertebra_labels.py').read())

Each label is a Draft Label: an annotation with a leader line and no shape,
drawn facing the camera. It is only for looking at. The Mesh Generator and the
Lead Designer read anatomy from STL files on disk (fem/scripts/config.py,
src/freecad/anatomy.py), never from the document, so the labels cannot reach a
mesh or a solve. They go in the "Labels" group: running again replaces the
vertebra labels and leaves anything else in that group alone, and Edit > Undo
takes a run back in one step. The document is not saved; save it from the GUI
if you want to keep them.

Vertebrae are found by their tree Label ("Vertebra T10"), which
body_aliases.yaml gives RADO's bodies and load_conformal.FCMacro carries over
to the extended set, so the same script works on either document.
"""
import re

import FreeCAD as App

GROUP = "SCS_Labels"
LABEL_RE = re.compile(r"^Vertebra (T\d+)$")
OFFSET_MM = 25.0    # text sits this far right (+x) of the transverse process tip
LEADER_MM = 10.0    # straight run of the leader next to the text
FONT = "20 mm"      # in Screen mode this is a fixed on-screen size
COLOUR = (1.0, 1.0, 1.0)


def vertebrae(doc):
    found = {}
    for obj in doc.Objects:
        m = LABEL_RE.match(obj.Label)
        if m and (hasattr(obj, "Mesh") or hasattr(obj, "Shape")):
            found[m.group(1)] = obj
    return found


def points(obj):
    if hasattr(obj, "Mesh"):
        return obj.Mesh.Points
    return [v.Point for v in obj.Shape.Vertexes]


def labels_group(doc, vertebra_objs):
    """The Labels group, emptied of earlier vertebra labels (those targeting a vertebra)."""
    group = doc.getObject(GROUP)
    if group is None:
        group = doc.addObject("App::DocumentObjectGroup", GROUP)
        group.Label = "Labels"
    for obj in list(group.Group):
        target = getattr(obj, "Target", None)
        if target and target[0] in vertebra_objs:
            doc.removeObject(obj.Name)
    return group


def main(doc=None):
    if not App.GuiUp:
        raise SystemExit("Draft labels are view objects; run this from the FreeCAD GUI")
    import Draft
    doc = doc or App.ActiveDocument
    found = vertebrae(doc)
    if not found:
        raise SystemExit("No 'Vertebra T<n>' bodies in %s" % doc.Label)
    doc.openTransaction("Vertebra labels")
    try:
        group = labels_group(doc, list(found.values()))
        for level, obj in sorted(found.items(), key=lambda item: int(item[0][1:])):
            tip = max(points(obj), key=lambda p: p.x)    # right transverse process
            tip = App.Vector(tip.x, tip.y, tip.z)
            label = Draft.make_label(target_point=tip, placement=tip + App.Vector(OFFSET_MM, 0, 0),
                                     target_object=obj, label_type="Custom", custom_text=level,
                                     direction="Horizontal", distance=-LEADER_MM)
            label.Label = "Label " + level
            view = label.ViewObject
            view.DisplayMode = "Screen"
            view.FontSize = FONT
            view.TextAlignment = "Middle"
            view.ArrowSizeStart = "1.5 mm"
            view.LineWidth = 2
            view.TextColor = COLOUR
            view.LineColor = COLOUR
            group.addObject(label)
        doc.recompute()
    except Exception:
        doc.abortTransaction()
        raise
    doc.commitTransaction()
    print("Labelled %s in %s" % (", ".join(sorted(found, key=lambda t: int(t[1:]))), doc.Label))


if __name__ == "__main__":
    main()
