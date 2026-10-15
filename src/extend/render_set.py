"""Render an anatomical STL set offscreen (pyvista), coloured by tissue.

    src/extend/.venv/bin/python src/extend/render_set.py DIR OUT.png [--tissues a,b] [--centre x,y,z]
        [--view x,y,z] [--zoom 1.0] [--opacity tissue=0.3,...] [--size 1400,1000]

Views without FreeCAD: useful remotely, and when FreeCAD is busy.
"""
import argparse
import os

import numpy as np

import meshkit as K

COLOR = {t["name"]: tuple(c / 255.0 for c in t["color_rgb"]) for t in K.ENTRIES if t.get("color_rgb")}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dir")
    ap.add_argument("png")
    ap.add_argument("--tissues", default=None, help="fem tissue names to show (default all)")
    ap.add_argument("--centre", default=None)
    ap.add_argument("--view", default="0,1,0", help="direction from the centre to the camera")
    ap.add_argument("--zoom", type=float, default=1.0)
    ap.add_argument("--opacity", default="", help="e.g. vertebra=0.25,epidural=0.15")
    ap.add_argument("--size", default="1400,1000")
    ap.add_argument("--title", default="")
    a = ap.parse_args()
    import pyvista as pv
    show = set(a.tissues.split(",")) if a.tissues else None
    opac = dict((kv.split("=")[0], float(kv.split("=")[1])) for kv in a.opacity.split(",") if kv)
    cls_of_pipe = {pipe: cls for pipe, cls in K.C.ANATOMY_CLASS}
    p = pv.Plotter(off_screen=True, window_size=[int(x) for x in a.size.split(",")])
    p.set_background("white")
    pts = []
    for fn in sorted(os.listdir(a.dir)):
        if not fn.lower().endswith(".stl"):
            continue
        pipe = K.pipe_of(os.path.splitext(fn)[0])
        if not pipe or (show and pipe not in show):
            continue
        V, T = K.weld(K.read_tri(os.path.join(a.dir, fn)), 0)
        m = pv.PolyData(V, np.hstack([np.full((len(T), 1), 3), T]).ravel())
        p.add_mesh(m, color=COLOR.get(cls_of_pipe[pipe], (0.7, 0.7, 0.7)), opacity=opac.get(pipe, 1.0),
                   smooth_shading=True, specular=0.3)
        pts.append(V.mean(0))
    c = np.array([float(x) for x in a.centre.split(",")]) if a.centre else np.mean(pts, 0)
    d = np.array([float(x) for x in a.view.split(",")])
    d /= np.linalg.norm(d)
    up = (0, 0, 1) if abs(d[2]) < 0.95 else (0, 1, 0)
    p.camera_position = [tuple(c + 300 * d), tuple(c), up]
    p.reset_camera_clipping_range()
    p.camera.zoom(a.zoom)
    if a.title:
        p.add_text(a.title, font_size=10, color="black")
    p.screenshot(a.png)
    p.close()
    print("wrote", a.png)


if __name__ == "__main__":
    main()
