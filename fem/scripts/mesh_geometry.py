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


def grid_points(lo, n, grid, start, stop):
    """Generate a slice of the structured grid without full coordinate arrays."""
    ijk = np.column_stack(np.unravel_index(np.arange(start, stop), tuple(n)))
    return np.asarray(lo) + grid * ijk
