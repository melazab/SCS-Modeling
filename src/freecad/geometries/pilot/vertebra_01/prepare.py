"""Measure one source vertebra and prepare a fitted profile and posterior mesh.

Run with Python containing numpy, scipy, trimesh, networkx, shapely and
mapbox_earcut. Source STL is opened read-only. Units are millimetres.
"""
from pathlib import Path
import hashlib
import json
import time

import numpy as np
from scipy.interpolate import splprep, splev
from scipy.spatial import cKDTree
import trimesh

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
SOURCE = ROOT / "STL_files" / "T8-10 - V1-2.STL"


def main():
    started = time.perf_counter()
    original = trimesh.load_mesh(SOURCE)
    mesh = original.copy()
    edges, counts = np.unique(mesh.edges_sorted, axis=0, return_counts=True)
    boundary = np.unique(edges[counts == 1])
    pairs = boundary[cKDTree(mesh.vertices[boundary]).query_pairs(2e-5, output_type="ndarray")]
    displacement = []
    for keep, discard in pairs:
        displacement.append(float(np.linalg.norm(mesh.vertices[keep] - mesh.vertices[discard])))
        mesh.faces[mesh.faces == discard] = keep
    mesh.remove_unreferenced_vertices()
    assert mesh.is_watertight and mesh.is_winding_consistent
    mesh.export(HERE / "seam_closed_reference.stl")

    normal = mesh.facets_normal[np.argmax(mesh.facets_area)].copy()
    normal[0] = 0
    normal /= np.linalg.norm(normal)
    if normal[2] < 0:
        normal *= -1
    rotation = np.array([[1, 0, 0], [0, normal[2], -normal[1]], normal])
    local = mesh.copy()
    local.vertices = mesh.vertices @ rotation.T
    cap_ids = np.argsort(mesh.facets_area)[-2:]
    cap_z = sorted(float(np.median(local.vertices[np.unique(mesh.faces[mesh.facets[i]]), 2])) for i in cap_ids)
    bottom, top = cap_z
    section = local.section(plane_normal=[0, 0, 1], plane_origin=[0, 0, bottom + 2])
    loops = [p[:, :2] for p in section.discrete if p[:, 1].mean() < 40 and np.ptp(p[:, 0]) > 20]
    assert len(loops) == 1
    outline = loops[0]
    arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(outline, axis=0), axis=1))]
    sample_arc = np.linspace(0, arc[-1], 513)
    samples = np.column_stack([np.interp(sample_arc, arc, outline[:, a]) for a in range(2)])
    # Fit the sectional polygon; the fit tolerance is not a whole-body error claim.
    tck, _ = splprep(samples.T, s=len(samples) * 0.015**2, per=True, k=3)
    knots, coeffs, degree = tck
    dense = np.array(splev(np.linspace(0, 1, 2001), tck)).T
    sampled_fit_error = cKDTree(dense).query(samples)[0]
    # Two approximately circular through-channels, fitted independently.
    channel_loops = [p[:, :2] for p in section.discrete if .9 < np.ptp(p[:, 0]) < 1.1]
    channels = []
    for p in channel_loops:
        c = np.linalg.lstsq(np.column_stack([2*p[:, 0], 2*p[:, 1], np.ones(len(p))]), (p*p).sum(1), rcond=None)[0]
        channels.append({"center": c[:2].tolist(), "radius": float(np.sqrt(c[2] + c[0]**2 + c[1]**2))})
    assert len(channels) == 2
    # Keep a narrow body overlap at the pedicle attachments; the final union is
    # measured separately. This split is geometric, not an anatomical boundary.
    split_y = 39.5
    posterior = local.slice_plane(plane_origin=[0, split_y, 0], plane_normal=[0, 1, 0], cap=True)
    assert posterior.is_watertight and posterior.is_winding_consistent
    posterior.vertices = posterior.vertices @ rotation
    posterior.export(HERE / "posterior_reference.stl")

    data = {
        "source_filename": SOURCE.name,
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_facets": len(original.faces), "source_vertices": len(original.vertices),
        "source_raw_watertight": bool(original.is_watertight),
        "source_boundary_edges": int((counts == 1).sum()),
        "source_volume_signed_mm3": float(original.volume),
        "source_bounds_mm": original.bounds.tolist(),
        "seam_weld_displacements_mm": displacement,
        "closed_reference_volume_mm3": float(mesh.volume),
        "closed_reference_euler": int(mesh.euler_number),
        "rotation_global_to_local": rotation.tolist(),
        "bottom_local_z_mm": bottom, "top_local_z_mm": top,
        "height_mm": top-bottom, "rim_radius_trial_mm": 1.5,
        "posterior_split_local_y_mm": split_y,
        "profile_samples_xy_mm": samples[:-1].tolist(),
        "spline_degree": degree, "spline_knots": knots.tolist(),
        "spline_coefficients_xy": np.array(coeffs).T.tolist(),
        "profile_sample_to_dense_spline_max_mm": float(sampled_fit_error.max()),
        "channels": sorted(channels, key=lambda p: p["center"][0]),
        "posterior_volume_mm3": float(posterior.volume),
        "preparation_wall_seconds": time.perf_counter()-started,
    }
    (HERE / "measurements.json").write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps({k:v for k,v in data.items() if k not in {"profile_samples_xy_mm", "spline_knots", "spline_coefficients_xy"}}, indent=2))


if __name__ == "__main__":
    main()
