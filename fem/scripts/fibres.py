"""White-matter fibre direction and the anisotropic conductivity tensor.

Khadka et al. 2020 give white matter 0.083 S/m across the fibres and 0.6 S/m
along them (config.SIGMA_ANISO, read from tissue_map.yaml). "Along" is the
cord's own axis, and the cord is curved: in the extended model its centreline
tangent turns from 20 deg one way at T12 to 23 deg the other way at T6. So the
direction is the local tangent of config.FIBRE_CENTRELINE at the nearest
centreline sample, not global z.

THE TENSOR. sigma = s_t I + a a^T, with a = sqrt(s_l - s_t) d and d the unit
fibre direction. Then sigma d = s_l d and sigma v = s_t v for v normal to d.

WHY BINS. Elmer takes one conductivity per material, so elements are grouped by
direction: consecutive centreline samples share a bin while their tangents stay
within BIN_DEG of the bin's first tangent, and each bin's direction is the
normalised mean of its members. Every element therefore carries a direction at
most BIN_DEG off its own nearest tangent; at 0.25 deg the longitudinal
conductivity leaks into the transverse direction by (s_l - s_t) sin^2 =
5e-6 S/m. The independent residual check applies the SAME binned tensors, so
it checks exactly what Elmer solved.

Output: (fibre_bin, fibre_a), with fibre_bin an int16 per element (-1 =
isotropic) and fibre_a the (K, 3) table of a-vectors.
"""
import json

import numpy as np

import config as C

BIN_DEG = 0.25
# White matter is within ~5 mm of the cord axis. A centroid farther than this
# from every centreline sample means the centreline does not cover the cord
# being meshed, and the direction would be a guess.
MAX_AXIS_DISTANCE_MM = 10.0


def load_centreline(path=None):
    """(points mm (n,3), unit tangents (n,3)) from extend_rado.py's samples."""
    with open(path or C.FIBRE_CENTRELINE) as fh:
        c = json.load(fh)
    p = np.asarray(c["points_mm"], float)
    t = np.asarray(c["tangents"], float)
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    return p, t


def direction_bins(tangents, bin_deg=BIN_DEG):
    """(sample_bin (n,), directions (K,3)): consecutive samples grouped while
    within bin_deg of the group's first tangent."""
    cos_max = np.cos(np.radians(bin_deg))
    sample_bin = np.empty(len(tangents), np.int64)
    dirs, start = [], 0
    for i in range(len(tangents)):
        if i > start and tangents[i] @ tangents[start] < cos_max:
            dirs.append(tangents[start:i].mean(0))
            start = i
        sample_bin[i] = len(dirs)
    dirs.append(tangents[start:].mean(0))
    dirs = np.array(dirs)
    return sample_bin, dirs / np.linalg.norm(dirs, axis=1, keepdims=True)


def tensor(transverse, a):
    """The 3x3 conductivity for one bin."""
    return transverse * np.eye(3) + np.outer(a, a)


def fibre_field(nodes, tets, lab, order, centreline=None, batch=500_000):
    """(fibre_bin int16 per tet, fibre_a (K,3), info) for the anisotropic
    classes in config.SIGMA_ANISO; (None, None, info) when there are none.

    `order` is the list `lab` indexes into. Bins are numbered per class, in
    ANISO order, so two anisotropic classes never share a row."""
    from scipy.spatial import cKDTree
    info = dict(anisotropic=dict(C.SIGMA_ANISO), bin_deg=BIN_DEG)
    present = [name for name in C.SIGMA_ANISO if name in order and np.any(lab == order.index(name))]
    if not present:
        return None, None, info
    points, tangents = load_centreline(centreline)
    sample_bin, dirs = direction_bins(tangents)
    tree = cKDTree(points)
    fibre_bin = np.full(len(tets), -1, np.int16)
    rows, worst = [], 0.0
    for name in present:
        transverse, longitudinal = C.SIGMA_ANISO[name]
        if longitudinal < transverse:
            raise ValueError("%s: longitudinal %.4g < transverse %.4g" % (name, longitudinal, transverse))
        offset = len(rows)
        rows.extend(np.sqrt(longitudinal - transverse) * dirs)
        sel = np.flatnonzero(lab == order.index(name))
        for start in range(0, len(sel), batch):
            part = sel[start:start + batch]
            centroid = nodes[tets[part]].mean(axis=1)
            dist, nearest = tree.query(centroid)
            worst = max(worst, float(dist.max()))
            fibre_bin[part] = offset + sample_bin[nearest]
    if len(rows) > np.iinfo(np.int16).max:
        raise ValueError("too many fibre bins (%d)" % len(rows))
    if worst > MAX_AXIS_DISTANCE_MM:
        raise ValueError("an anisotropic element is %.1f mm from the cord centreline (%s); "
                         "the centreline does not cover this cord" % (worst, C.FIBRE_CENTRELINE))
    used = np.unique(fibre_bin[fibre_bin >= 0])
    info.update(classes=present, bins=int(len(rows)), bins_used=int(len(used)),
                max_axis_distance_mm=round(worst, 3),
                turn_deg=round(float(np.degrees(np.arccos(np.clip(dirs[0] @ dirs[-1], -1, 1)))), 2))
    return fibre_bin, np.asarray(rows), info
