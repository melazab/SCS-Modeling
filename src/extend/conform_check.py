"""Check a conform.py output set against its source set.

    src/extend/.venv/bin/python src/extend/conform_check.py [--src DIR] [--out DIR] [--n 400000]

Random points, half spread over the model and half within 0.3 mm of a source
surface (where interfaces are), are labelled two ways:

  * source set, the way fem/ labels a tet: the first tissue in the FEM's
    priority order with any body containing the point;
  * output set: the tissue of the output body containing it. More than one
    containing body is an overlap and is counted.

Agreement should be ~100% everywhere except inside the blend fillets, which
add a little material where a joint's crease was.
"""
import argparse
import json
import os
import time

import numpy as np

import meshkit as K


def bodies_of(d, skip_trunk_segments=False):
    out = []
    for fn in sorted(os.listdir(d)):
        if fn.lower().endswith(".stl"):
            pipe = K.pipe_of(os.path.splitext(fn)[0])
            if pipe:
                V, T = K.stlio.read_stl(os.path.join(d, fn))
                out.append((fn, pipe, K.InsideTester(V, T), V.min(0), V.max(0)))
    return out


def containment(pts, bodies, chunk=20000):
    """Per point: list of indices of bodies containing it (points sorted in z
    so each chunk's box is small and most bodies are culled)."""
    order = np.argsort(pts[:, 2])
    hits = [[] for _ in range(len(pts))]
    for s in range(0, len(pts), chunk):
        idx = order[s:s + chunk]
        P = pts[idx]
        lo, hi = P.min(0), P.max(0)
        for bi, (_, _, t, blo, bhi) in enumerate(bodies):
            if np.any(bhi < lo) or np.any(blo > hi):
                continue
            m = np.all((P >= blo) & (P <= bhi), 1)
            if not m.any():
                continue
            ins = np.zeros(len(P), bool)
            ins[m] = t(P[m])
            for i in idx[ins]:
                hits[i].append(bi)
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", default=os.path.join(K.C.OUT, "geometry", "rado_up4"))
    ap.add_argument("--out", default=os.path.join(K.C.OUT, "geometry", "rado_up4_conformal"))
    ap.add_argument("--n", type=int, default=400000)
    a = ap.parse_args()
    t0 = time.time()
    src, dst = bodies_of(a.src), bodies_of(a.out)
    rng = np.random.default_rng(0)
    lo = np.min([b[3] for b in dst], 0)
    hi = np.max([b[4] for b in dst], 0)
    uni = lo + rng.random((a.n // 2, 3)) * (hi - lo)
    # near-surface points: jittered vertices of the source bodies inside the output's box
    allv = np.vstack([K.stlio.read_stl(os.path.join(a.src, b[0]))[0] for b in src])
    allv = allv[np.all((allv >= lo) & (allv <= hi), 1)]
    near = allv[rng.integers(0, len(allv), a.n - a.n // 2)] + rng.normal(scale=0.15, size=(a.n - a.n // 2, 3))
    pts = np.vstack([uni, near])
    order = K.ORDER
    h_src = containment(pts, src)
    h_dst = containment(pts, dst)
    lab_src = np.array([min((order.index(src[b][1]) for b in h), default=-1) for h in h_src])
    n_dst = np.array([len(h) for h in h_dst])
    lab_dst = np.array([order.index(dst[h[0]][1]) if len(h) else -1 for h in h_dst])
    res = {"points": len(pts), "uniform": len(uni), "near_surface": len(near),
           "points_in_more_than_one_output_body": int((n_dst > 1).sum())}
    for name, sl in (("all", slice(None)), ("uniform", slice(0, len(uni))), ("near_surface", slice(len(uni), None))):
        agree = lab_src[sl] == lab_dst[sl]
        res["agreement_" + name] = round(float(agree.mean()), 6)
    dis = lab_src != lab_dst
    pairs = {}
    for s_, d_ in zip(lab_src[dis], lab_dst[dis]):
        key = "%s -> %s" % (order[s_] if s_ >= 0 else "background", order[d_] if d_ >= 0 else "background")
        pairs[key] = pairs.get(key, 0) + 1
    res["disagreements_source_to_output"] = dict(sorted(pairs.items(), key=lambda kv: -kv[1]))
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(a.out, "conform_check.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
