"""Shared anatomical coverage and size regions for meshing and cost estimation."""
import numpy as np
import config as C
from stlio import read_stl

# These are mesh controls, not material assignments. Conductivities and
# overlapping-body precedence remain defined by config/tissue_map.
REGION_CLASSES = {
    'dura': ('dura',),
    'white': ('white', 'grey'),
    'canal': ('epidural', 'csf'),
    'fine': ('root', 'drg', 'blood', 'sympathetic'),
    'coarse': ('vertebra', 'disc'),
}
TARGET = dict(dura=0.45, white=0.70, canal=0.90, fine=0.45, coarse=2.0, lead=0.25)


def geometry(contact_stl=None, insulator_stl=None, factor=1.0):
    from build_mesh import densify
    contact_stl = C.CONTACT_STL if contact_stl is None else contact_stl
    insulator_stl = C.INSULATOR_STL if insulator_stl is None else insulator_stl
    clouds, bounds = {}, []
    for region, classes in REGION_CLASSES.items():
        pieces = []
        for cls in classes:
            for path in C.TISSUE_BODIES[cls]:
                v, t = read_stl(path)
                bounds.extend((v.min(0), v.max(0)))
                pieces.append(densify(v, t, TARGET[region] * factor))
        clouds[region] = np.concatenate(pieces)
    pieces = []
    for path in list(contact_stl.values()) + [insulator_stl]:
        v, t = read_stl(path)
        bounds.extend((v.min(0), v.max(0)))
        pieces.append(densify(v, t, TARGET['lead'] * factor))
    clouds['lead'] = np.concatenate(pieces)
    return clouds, np.asarray(bounds)


SLAB_MM = 1.0


class Extent(object):
    """Where the anatomy is, slab by slab along z, for crop windows.

    Built once from the size-field clouds (densified surfaces of every body
    and the lead), so a cropped box hugs the anatomy that is actually inside
    the z window instead of the whole model's x/y extent -- the extended
    model's canal moves ~25 mm in y between T12 and T6."""

    def __init__(self, clouds, bounds):
        self.lo_all, self.hi_all = bounds.min(0), bounds.max(0)
        self.slab_lo = None
        if not clouds or 'lead' not in clouds:
            return                    # whole-model boxes only
        pts = np.concatenate([c for c in clouds.values()])
        self.z0 = float(pts[:, 2].min())
        k = ((pts[:, 2] - self.z0) // SLAB_MM).astype(np.int64)
        n = int(k.max()) + 1
        self.slab_lo = np.full((n, 3), np.inf)
        self.slab_hi = np.full((n, 3), -np.inf)
        order = np.argsort(k, kind='stable')
        k, pts = k[order], pts[order]
        starts = np.searchsorted(k, np.arange(n))
        ends = np.searchsorted(k, np.arange(n), side='right')
        for i in np.flatnonzero(ends > starts):
            self.slab_lo[i] = pts[starts[i]:ends[i]].min(0)
            self.slab_hi[i] = pts[starts[i]:ends[i]].max(0)
        lead = clouds['lead']
        self.lead_z = (float(lead[:, 2].min()), float(lead[:, 2].max()))

    def box(self, crop_mm, margin):
        """(lo, hi) of the mesh box. crop_mm <= 0: all anatomy + margin, as
        always. Otherwise z runs from crop_mm below the lead's lowest point to
        crop_mm above its highest (never past the anatomy), x and y cover
        what lies in that window, plus margin. The z faces cut straight
        through tissue; they get no margin."""
        lo, hi = self.lo_all - margin, self.hi_all + margin
        if not crop_mm or crop_mm <= 0:
            return lo, hi
        if self.slab_lo is None:
            raise ValueError('a crop window needs the anatomy and lead surface clouds')
        zlo = max(self.lead_z[0] - crop_mm, self.lo_all[2])
        zhi = min(self.lead_z[1] + crop_mm, self.hi_all[2])
        i0 = max(0, int((zlo - self.z0) // SLAB_MM))
        i1 = min(len(self.slab_lo), int((zhi - self.z0) // SLAB_MM) + 1)
        wlo = self.slab_lo[i0:i1].min(0) - margin
        whi = self.slab_hi[i0:i1].max(0) + margin
        return (np.array([wlo[0], wlo[1], zlo]), np.array([whi[0], whi[1], zhi]))


def grid_points(lo, n, grid, start, stop):
    """Generate a slice of the structured grid without full coordinate arrays."""
    ijk = np.column_stack(np.unravel_index(np.arange(start, stop), tuple(n)))
    return np.asarray(lo) + grid * ijk
