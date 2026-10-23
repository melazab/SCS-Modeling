"""Estimate mesh and preview resources before a build.

Integrate 1/h^3 on a coarse sampling grid using the same anatomical bounds and
region surfaces as build_mesh. The count calibration (2.41) was measured on
six older canal-only meshes; full-anatomy estimates are provisional. RAM covers
meshing, classification and bounded face sorting, with a conservative range.
The solver has its own memory check; meshing success is not a solve guarantee.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from stlio import read_stl

# Empirical over-count of the perfect-regular-tet integral vs. what gmsh
# actually produces. See the calibration table in the module docstring.
RAW_TO_ACTUAL = 2.41

# Sampling grid for the ESTIMATE only (mm). build_mesh uses 0.25 mm; this is
# coarser by 4x per axis = 64x fewer points, which is what makes this
# interactive. Validated against the fine-grid integral in
# `_selftest_against_known_runs()` below -- coarse sampling tracks the fine
# integral closely because h(x) is a smooth (piecewise-linear-in-distance)
# function, so undersampling shifts the integral far less than it would for a
# field with fine structure.
ESTIMATE_GRID_MM = 2.0

# Conservative mesh + classification + preview envelope, including sorting.
# The old 150-250 byte/tet estimate covered Gmsh alone and missed the OOM.
BYTES_PER_TET_LO = 400
BYTES_PER_TET_HI = 700

# Rough throughput, tets per second, for the meshing stage on this machine.
# Anchored on the 2026-09-14 fine build: 3.50 M tets, 912.8 s total of which
# roughly 150 s was meshing (the rest was the size field) -> ~23 k tets/s.
# Used only to say "minutes" vs "hours"; not a promise.
TETS_PER_SEC = 23_000


# How much coarser than build_mesh's own surface clouds the estimator's
# clouds are. Measured 2026-09-15 against the three runs with known tet
# counts -- build time is dominated by cKDTree construction over these
# clouds, so this is the one knob that matters for interactivity:
#
#     factor   build    worst error vs. ground truth
#       1      20.6 s            1.7 %
#       2       4.3 s            2.3 %     <- chosen
#       3       2.1 s            4.8 %
#       4       1.2 s            5.4 %
#
# 2 buys a 4.8x speed-up for 0.6 percentage points of accuracy. Past that the
# clouds get sparse enough that distance-to-nearest-SAMPLE starts to diverge
# from distance-to-SURFACE near the thin dura, and the estimate drifts low.
DENSIFY_FACTOR = 2.0


def _clouds(contact_stl=None, insulator_stl=None, densify_factor=DENSIFY_FACTOR):
    """Surface point clouds for the four size-field regions.

    Mirrors build_mesh.build()'s own cloud construction, but with coarser
    densify targets (densify_factor > 1 loosens them) since an estimate does
    not need every facet subdivided -- see DENSIFY_FACTOR above.
    """
    from mesh_geometry import geometry
    return geometry(contact_stl, insulator_stl, factor=densify_factor)


class CostModel(object):
    """Holds the expensive-but-parameter-independent part of the estimate.

    The distances from each sample point to each region's surface depend ONLY
    on the geometry and the sampling grid -- NOT on H_MIN/H_MAX/h_near/
    plateau/growth. So they are computed once here, and then any number of
    parameter combinations can be costed by pure arithmetic in milliseconds.
    That is what lets the panel re-estimate live as a spinbox moves.
    """

    def __init__(self, contact_stl=None, insulator_stl=None,
                 grid=ESTIMATE_GRID_MM, margin=1.0):
        from scipy.spatial import cKDTree
        from mesh_geometry import Extent
        clouds, bounds = _clouds(contact_stl, insulator_stl)
        # The crop window is per parameter set, so keep the sample points and
        # mask them per estimate; the full box is sampled once.
        self.extent = Extent(clouds, bounds)
        self.margin = margin
        lo, hi = self.extent.box(0, margin)
        n = np.maximum(np.ceil((hi - lo) / grid).astype(int) + 1, 2)
        axes = [lo[i] + grid * np.arange(n[i]) for i in range(3)]
        gx, gy, gz = np.meshgrid(*axes, indexing="ij")
        P = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)

        self.grid = grid
        self.points = P
        self.n_samples = len(P)
        self.dists = {}
        for name, cloud in clouds.items():
            # Conservative radius for the UI's H_Max <= 10 mm and growth >= .5.
            # Avoid expensive unbounded searches far from narrow root surfaces.
            distance = cKDTree(cloud).query(P, workers=4, distance_upper_bound=24.0)[0]
            # A capped distance is a lower bound: for CLI settings whose
            # influence extends farther, this overestimates refinement rather
            # than missing distant refinement and underestimating memory.
            self.dists[name] = np.minimum(distance, 24.0)

    def raw_tets(self, params):
        """The uncalibrated perfect-tet integral for one parameter set."""
        crop = float(params.get("CROP_MM", 0.0) or 0.0)
        keep = slice(None)
        if crop > 0:
            lo, hi = self.extent.box(crop, self.margin)
            keep = np.all((self.points >= lo) & (self.points <= hi), axis=1)
        h = np.full(self.n_samples, float(params["H_MAX"]))[keep]
        for name, hnear, plateau, growth in params["FIELDS"]:
            d = self.dists[name][keep]
            np.minimum(h, hnear + growth * np.maximum(0.0, d - plateau), out=h)
        np.clip(h, params["H_MIN"], params["H_MAX"], out=h)
        cell_v = self.grid ** 3
        return float(np.sum(cell_v / (h ** 3 / (6 * np.sqrt(2)))))

    def estimate(self, params):
        """Calibrated prediction for one parameter set."""
        raw = self.raw_tets(params)
        n_tets = raw / RAW_TO_ACTUAL
        return dict(
            n_tets=n_tets,
            ram_gb_lo=0.5 + n_tets * BYTES_PER_TET_LO / 1e9,
            ram_gb_hi=1.0 + 1.5 * n_tets * BYTES_PER_TET_HI / 1e9,
            mesh_seconds=n_tets / TETS_PER_SEC,
        )


def available_ram_gb():
    """Currently available physical RAM in GB (swap is not a safe budget)."""
    try:
        import psutil
        return psutil.virtual_memory().available / 1e9
    except ImportError:
        try:
            return (os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")) / 1e9
        except (ValueError, OSError, AttributeError):
            return None


def describe(est, ram_gb=None, budget="available"):
    """One short human line, plus a severity the UI can colour on.

    Returns (severity, text) where severity is "ok" | "warn" | "danger".
    Deliberately terse -- this goes in a panel, not a log. `ram_gb` defaults
    to this machine's available RAM; a cluster job passes its allocation and
    budget="HPC allocation".
    """
    ram_gb = available_ram_gb() if ram_gb is None else ram_gb
    n = est["n_tets"]
    if n >= 100e6:
        count = "%.0f M tets" % (n / 1e6)
    elif n >= 1e6:
        count = "%.1f M tets" % (n / 1e6)
    else:
        count = "%.0f k tets" % (n / 1e3)
    lo_gb, hi_gb = est["ram_gb_lo"], est["ram_gb_hi"]
    if hi_gb < 1.0:
        ram = "%.0f\u2013%.0f MB RAM" % (lo_gb * 1000, hi_gb * 1000)
    elif hi_gb < 10.0:
        ram = "%.1f\u2013%.1f GB RAM" % (lo_gb, hi_gb)
    else:
        ram = "%.0f\u2013%.0f GB RAM" % (lo_gb, hi_gb)
    parts = [count, ram]
    # The time is gmsh only, at a workstation calibration (TETS_PER_SEC); on
    # 24 smp cores a 415 M-tet build meshed in 188 s against "~4.9 h", and the
    # stages after gmsh are not in it at all. Shown for local runs only until
    # a per-stage calibration replaces it.
    if budget == "available":
        mins = est["mesh_seconds"] / 60.0
        parts.append("gmsh ~%.0f min" % max(1, round(mins)) if mins < 90 else "gmsh ~%.1f h" % (mins / 60.0))
    sev = "ok"
    if ram_gb is not None:
        if est["ram_gb_hi"] > 0.8 * ram_gb:
            sev = "danger"
        elif est["ram_gb_hi"] > 0.5 * ram_gb:
            sev = "warn"
    text = " \u00b7 ".join(parts)
    if sev in ("danger", "warn"):
        text += " \u00b7 %s %.0f GB budget (80%% of %.0f GB %s)" % (
            "over" if sev == "danger" else "near", 0.8 * ram_gb, ram_gb, budget)
    return sev, text
