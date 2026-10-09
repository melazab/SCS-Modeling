"""Bounded-memory P1 element calculations, assembly and residual checks."""
import os
import numpy as np
import scipy.sparse as sp


def batch_size():
    # Approximately 1 KiB of element temporaries; one batch at a time.
    memory = float(os.environ.get('SCS_SOLVE_MEMORY_GB', '4'))
    return max(1, min(100_000, int(memory * 1e9 * .02 / 1024)))


def elements(nodes, tets, sigma, batch=None):
    """Yield (start, positive volumes in m^3, local stiffness matrices).

    Signed Jacobian inversion handles either orientation. Absolute volume
    supplies the integration weight; connectivity need not be copied/flipped.
    """
    batch = batch_size() if batch is None else int(batch)
    if batch < 1:
        raise ValueError('Element batch size must be positive')
    for start in range(0, len(tets), batch):
        tet = tets[start:start + batch]
        p = nodes[tet] * 1e-3
        e1, e2, e3 = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]
        # Closed-form inverse of the Jacobian [e1 e2 e3]: its rows are the
        # cross products over the determinant. ~3x faster than np.linalg's
        # per-matrix LAPACK calls, which dominated the residual check.
        g = np.empty((len(tet), 4, 3))
        g[:, 1], g[:, 2], g[:, 3] = np.cross(e2, e3), np.cross(e3, e1), np.cross(e1, e2)
        det = np.einsum('ij,ij->i', e1, g[:, 1])
        if not np.isfinite(det).all() or np.any(det == 0):
            raise ValueError('Nonfinite or degenerate tetrahedron in element batch')
        volume = np.abs(det) / 6
        g[:, 1:] /= det[:, None, None]
        g[:, 0] = -g[:, 1:].sum(axis=1)
        ke = np.einsum('eik,ejk->eij', g, g)
        ke *= (sigma[start:start + len(tet)] * volume)[:, None, None]
        yield start, volume, ke


def assemble(nodes, tets, sigma, batch=None):
    """Reference sparse assembly without global element matrices/triplets."""
    matrix = sp.csr_matrix((len(nodes), len(nodes)), dtype=np.float64)
    volumes = np.empty(len(tets))
    for start, vol, ke in elements(nodes, tets, sigma, batch):
        tet = tets[start:start + len(vol)]
        volumes[start:start + len(vol)] = vol
        rows = np.repeat(tet, 4, axis=1).ravel()
        cols = np.tile(tet, (1, 4)).ravel()
        part = sp.coo_matrix((ke.ravel(), (rows, cols)), shape=matrix.shape).tocsr()
        matrix = matrix + part
    return matrix, volumes, tets


def action(nodes, tets, sigma, values, batch=None):
    """Compute A @ values without constructing the global sparse matrix."""
    result = np.zeros(len(nodes))
    for start, vol, ke in elements(nodes, tets, sigma, batch):
        tet = tets[start:start + len(vol)]
        local = np.einsum('eij,ej->ei', ke, values[tet])
        np.add.at(result, tet.ravel(), local.ravel())
    return result
