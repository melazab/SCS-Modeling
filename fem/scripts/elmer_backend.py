"""Elmer StatCurrentSolver backend with bounded preparation and verification.

One unit-current basis per actual lead contact, uniform volumetric injection,
uniform-area boundary return, and one pinned gauge. Result schema matches the
visualizer. No Python global stiffness matrix is needed in production.
"""
import glob
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import config as C
import fibres
from element_batches import elements, action, batch_size
from mesh_faces import face_neighbours, face_chunks

# SCS_ELMER_BIN points at another install, e.g. the Pioneer build in ~/opt.
ELMER_BIN = Path(os.environ.get('SCS_ELMER_BIN', '/opt/elmerfem/bin'))


def run(command, cwd, log):
    env = dict(os.environ)
    env['ELMER_HOME'] = str(ELMER_BIN.parent)
    env['ELMER_SOLVER_HOME'] = str(ELMER_BIN.parent / 'share/elmersolver')
    dependencies = ELMER_BIN.parent / 'deps/usr/lib/x86_64-linux-gnu'
    if dependencies.is_dir():
        env['LD_LIBRARY_PATH'] = str(dependencies) + os.pathsep + env.get('LD_LIBRARY_PATH', '')
    # CPUs are MPI tasks. Do not multiply this allocation by BLAS/OpenMP threads.
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        env[name] = '1'
    print('Running: ' + ' '.join(map(str, command)), flush=True)
    with open(log, 'w') as output:
        status = subprocess.run(list(map(str, command)), cwd=cwd, env=env,
                                stdout=output, stderr=subprocess.STDOUT)
    if status.returncode:
        with open(log, 'rb') as output:
            output.seek(max(0, os.path.getsize(log) - 5000))
            tail = output.read().decode(errors='replace')
        raise RuntimeError('Elmer command failed (exit %d). Log: %s\n%s' %
                           (status.returncode, log, tail))


def connectivity(tets, n):
    """Reject disconnected active domains (tets joined through shared nodes).

    One sparse build; summing per-batch CSR matrices was quadratic in the
    number of batches (~520 at 52 M tets)."""
    row = np.repeat(tets[:, 0], 3)
    col = tets[:, 1:].ravel()
    graph = sp.csr_matrix((np.ones(len(row), bool), (row, col)), shape=(n, n))
    del row, col
    count = connected_components(graph, directed=False, return_labels=False)
    if count != 1:
        raise RuntimeError('Active mesh has %d disconnected components; check tissue connectivity before solving.' % count)


def element_tags(labels, order, fibre_bin=None):
    """Elmer body tag per element, numbered 1..N over the bodies present.

    One body per tissue label, one for background, and one per anisotropic
    fibre bin, so every bin carries its own conductivity tensor. Elmer
    requires Body 1..N with no gap ("Entry missing for: Body 3"), so the raw
    keys (label, background = len(order), bin = len(order) + 1 + bin) are
    renumbered over the ones that actually occur. A tissue with no elements --
    white matter whose elements all went to fibre bins, or a class a crop
    window leaves out -- would otherwise leave a hole."""
    key = np.where(labels >= 0, labels, len(order)).astype(np.int32)
    if fibre_bin is not None:
        m = fibre_bin >= 0
        key[m] = len(order) + 1 + fibre_bin[m]
    _, tag = np.unique(key, return_inverse=True)
    return (tag + 1).astype(np.int32)


# Pieces of the active mesh that are not connected to the lead are dropped
# before solving, up to this fraction of the active volume. They carry no
# current from the lead, and as separate pieces they make the system singular.
# On the fused T6-T12 set they are specks of peripheral vessel tips, isolated
# where a vessel is thinner than the elements (12 pieces, 9 mm3 of 462 cm3 at
# coarse settings, 2026-10-15). More than this means a real modelling gap, so
# the solve refuses instead.
ISLAND_MAX_FRACTION = 0.005


def _volumes(nodes, tets, chunk=2_000_000):
    out = np.empty(len(tets))
    for start in range(0, len(tets), chunk):
        p = nodes[tets[start:start + chunk]]
        out[start:start + len(p)] = np.abs(np.einsum(
            'ij,ij->i', p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))) / 6
    return out


def drop_islands(nodes, tets, labels, order):
    """(keep mask over tets, report): the connected piece(s) holding the lead's
    contacts are kept, everything else dropped and reported by tissue.
    Connectivity is through shared nodes, as in connectivity()."""
    used, local = np.unique(tets, return_inverse=True)
    local = local.reshape(tets.shape)
    n = len(used)
    row = np.repeat(local[:, 0], 3)
    col = local[:, 1:].ravel()
    graph = sp.csr_matrix((np.ones(len(row), bool), (row, col)), shape=(n, n))
    del row, col
    count, comp = connected_components(graph, directed=False)
    del graph
    report = dict(components=int(count), kept_components=int(count), dropped_tets=0,
                  dropped_mm3=0.0, dropped_by_tissue={})
    if count == 1:
        return np.ones(len(tets), bool), report
    tet_comp = comp[local[:, 0]]
    del local, comp
    contacts = [order.index(name) for name in order if name.startswith('contact')]
    lead = np.unique(tet_comp[np.isin(labels, contacts)])
    keep = np.isin(tet_comp, lead)
    vol = _volumes(nodes, tets)
    dropped_lab = labels[~keep]
    dropped_vol = vol[~keep]
    report.update(kept_components=int(len(lead)), dropped_tets=int((~keep).sum()),
                  dropped_mm3=float(dropped_vol.sum()),
                  dropped_by_tissue={(order[k] if k >= 0 else 'background'):
                                     float('%.4g' % dropped_vol[dropped_lab == k].sum())
                                     for k in np.unique(dropped_lab)})
    if dropped_vol.sum() > ISLAND_MAX_FRACTION * vol.sum():
        raise RuntimeError('%d mesh pieces (%.1f mm3, %s) are not connected to the lead -- more than %.1f%% '
                           'of the domain; check tissue connectivity before solving.'
                           % (count - len(lead), dropped_vol.sum(), report['dropped_by_tissue'],
                              100 * ISLAND_MAX_FRACTION))
    return keep, report


def write_mesh(directory, nodes, tets, labels, order, pin=None, nbr=None, fibre_bin=None, tags=None):
    """Write native Elmer mesh in SI units; return normalized nodal return."""
    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / 'mesh.nodes', 'w') as stream:
        for start in range(0, len(nodes), batch_size()):
            xyz = nodes[start:start + batch_size()] * 1e-3
            ids = np.arange(start + 1, start + len(xyz) + 1)
            np.savetxt(stream, np.column_stack((ids, np.full(len(ids), -1), xyz)),
                       fmt=['%d', '%d', '%.17g', '%.17g', '%.17g'])
    tags = element_tags(labels, order, fibre_bin) if tags is None else tags
    with open(directory / 'mesh.elements', 'w') as stream:
        for start in range(0, len(tets), batch_size()):
            tet = tets[start:start + batch_size()]
            tag = tags[start:start + len(tet)]
            ids = np.arange(start + 1, start + len(tet) + 1)
            np.savetxt(stream, np.column_stack((ids, tag, np.full(len(tet), 504), tet + 1)), fmt='%d')
    nodal_area = np.zeros(len(nodes))
    if pin is None:
        pin = int(np.argmax(np.linalg.norm(nodes[:, :2] - np.array([55.60, 87.0]), axis=1)))
    count = 0
    with open(directory / 'mesh.boundary', 'w') as stream:
        nbr = face_neighbours(tets) if nbr is None else nbr
        for faces, owner, other in face_chunks(tets, nbr):
            lone = other < 0
            boundary = faces[lone]
            p = nodes[boundary] * 1e-3
            area = np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1) / 2
            np.add.at(nodal_area, boundary.ravel(), np.repeat(area / 3, 3))
            ids = np.arange(count + 1, count + len(boundary) + 1)
            records = np.column_stack((ids, np.ones(len(ids), int), owner[lone] + 1,
                                       np.zeros(len(ids), int), np.full(len(ids), 303), boundary + 1))
            np.savetxt(stream, records, fmt='%d')
            count += len(boundary)
        # Point boundary elements are correctly remapped by ElmerGrid/MPI.
        # Target Nodes in a SIF is local numbering and would pin one WRONG
        # node per rank after partitioning.
        parent = None
        for start in range(0, len(tets), batch_size()):
            hits = np.flatnonzero(np.any(tets[start:start + batch_size()] == pin, axis=1))
            if len(hits):
                parent = start + int(hits[0]) + 1
                break
        if parent is None:
            raise ValueError('Gauge node is not part of the mesh')
        stream.write('%d 2 %d 0 101 %d\n' % (count + 1, parent, pin + 1))
    (directory / 'mesh.header').write_text('%d %d %d\n3\n504 %d\n303 %d\n101 1\n' %
                                         (len(nodes), len(tets), count + 1, len(tets), count))
    total = float(nodal_area.sum())
    if total <= 0:
        raise ValueError('Mesh has no exterior boundary area')
    return nodal_area / total, total


def sif(order, labels, sigma, contact, source_volume, return_area, pin, fibre=None, tags=None):
    if tags is None:
        tags = element_tags(labels, order, None if fibre is None else fibre[0])
    bodies = []
    for tag, first in zip(*np.unique(tags, return_index=True)):
        label = int(labels[first])
        name = order[label] if label >= 0 else 'background'
        if fibre is not None and fibre[0][first] >= 0:
            tensor = fibres.tensor(float(sigma[first]), fibre[1][fibre[0][first]])
            cond = 'Electric Conductivity(3,3) = Real ' + ' '.join('%.17g' % v for v in tensor.ravel())
        else:
            cond = 'Electric Conductivity = %.17g' % float(sigma[first])
        bodyforce = '\n  Body Force = 1' if name == 'contact%d' % contact else ''
        bodies.append('Body %d\n  Target Bodies(1) = %d\n  Equation = 1\n  Material = %d%s\nEnd\n'
                      'Material %d\n  %s\nEnd\n' %
                      (tag, tag, tag, bodyforce, tag, cond))
    return '''Header
  CHECK KEYWORDS Warn
  Mesh DB "." "mesh"
  Results Directory "results"
End
Simulation
  Max Output Level = 5
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Post File = "field.vtu"
End
Equation 1
  Active Solvers(1) = 1
End
Solver 1
  Equation = Stat Current Solver
  Procedure = "StatCurrentSolve" "StatCurrentSolver"
  Variable = Potential
  Variable DOFs = 1
  Calculate Volume Current = False
  Calculate Joule Heating = False
  Linear System Solver = Iterative
  Linear System Use Hypre = True
  Linear System Iterative Method = CG
  Linear System Symmetric = True
  Linear System Preconditioning = BoomerAMG
  BoomerAMG Coarsen Type = 10
  BoomerAMG Interpolation Type = 6
  BoomerAMG Relax Type = 6
  BoomerAMG Strong Threshold = 0.5
  Linear System Max Iterations = 4000
  Linear System Convergence Tolerance = 1.0e-11
  Linear System Abort Not Converged = True
  Linear System Residual Output = 100
  Steady State Convergence Tolerance = 1.0e-8
End
Body Force 1
  Current Source = %.17g
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Current Density BC = True
  Current Density = %.17g
End
Boundary Condition 2
  Target Boundaries(1) = 2
  Potential = 0.0
End
%s
''' % (1. / source_volume, -1. / return_area, '\n'.join(bodies))


def read_potential(results, nodes, tree=None, workers=1):
    """Map partition output to original nodes and check interface agreement.

    Pass `tree` (cKDTree of nodes in metres) to reuse it across contacts."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    tree = cKDTree(nodes * 1e-3) if tree is None else tree
    values = np.full(len(nodes), np.nan)
    files = sorted(results.glob('*.vtu'))
    if not files:
        raise RuntimeError('Elmer produced no VTU results')
    for path in files:
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(path)); reader.Update()
        grid = reader.GetOutput()
        arr = grid.GetPointData().GetArray('potential')
        if arr is None:
            arr = grid.GetPointData().GetArray('Potential')
        if arr is None or grid.GetNumberOfPoints() == 0:
            raise RuntimeError('Missing potential data in %s' % path)
        v = vtk_to_numpy(arr)
        dist, ids = tree.query(vtk_to_numpy(grid.GetPoints().GetData()), workers=workers)
        if not np.isfinite(v).all() or np.max(dist) > 1e-9:
            raise RuntimeError('Elmer result coordinates/values do not match input mesh')
        present = np.isfinite(values[ids])
        if present.any() and not np.allclose(values[ids[present]], v[present], rtol=1e-8, atol=1e-8):
            raise RuntimeError('MPI partitions disagree at shared nodes')
        values[ids] = v
    if not np.isfinite(values).all():
        raise RuntimeError('Elmer result is missing nodes')
    return values


def accepted(rel, reaction):
    from assign_and_solve import ACCEPTED_RESIDUAL
    return bool(np.isfinite(rel) and rel <= ACCEPTED_RESIDUAL and abs(reaction) <= 1e-6)


def finished(case, text):
    """An earlier Elmer run of exactly this case.sif ran to the end."""
    try:
        return ((case / 'case.sif').read_text() == text
                and 'ALL DONE' in (case / 'solver.log').read_text(errors='replace')
                and any((case / 'results').glob('*.vtu')))
    except OSError:
        return False


def solve_prepared(nodes, tets, labels, sigma, order, contact_ids, directory, resources,
                   progress=lambda p: None, nbr=None, fibre=None):
    """Solve a classified mesh; independently verify every unit-current field.

    nbr is the face adjacency of tets (mesh_faces.face_neighbours), if the
    caller already has it. fibre is (fibre_bin, fibre_a) from
    fibres.fibre_field() for anisotropic tissue, or None."""
    from assign_and_solve import ACCEPTED_RESIDUAL
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name in ('ElmerSolver_mpi', 'ElmerGrid'):
        if not (ELMER_BIN / name).is_file():
            raise RuntimeError('Required Elmer executable missing: %s' % (ELMER_BIN / name))
    mpi = shutil.which('mpirun')
    if resources.cpus > 1 and not mpi:
        raise RuntimeError('mpirun not found for local parallel solve')
    resources.mark('connectivity')
    connectivity(tets, len(nodes))
    resources.mark('element preparation')
    vol = np.empty(len(tets))
    for start, v, _ in elements(nodes, tets, sigma):
        vol[start:start + len(v)] = v
    pin = int(np.argmax(np.linalg.norm(nodes[:, :2] - np.array([55.60, 87.0]), axis=1)))
    resources.mark('mesh export')
    if nbr is None:
        nbr = face_neighbours(tets, resources.cpus)
    tags = element_tags(labels, order, None if fibre is None else fibre[0])
    ret, area = write_mesh(directory / 'mesh', nodes, tets, labels, order, pin=pin, nbr=nbr, tags=tags)
    del nbr
    progress(12)
    if resources.cpus > 1:
        resources.mark('mesh partitioning')
        run([ELMER_BIN / 'ElmerGrid', '2', '2', 'mesh', '-metis', str(resources.cpus), '-parthypre'],
            directory, directory / 'partition.log')
    phi = np.empty((len(contact_ids), len(nodes)))
    residuals = []
    reactions = []
    free = np.ones(len(nodes), bool); free[pin] = False
    version = tree = None
    for pos, contact in enumerate(contact_ids):
        mask = labels == order.index('contact%d' % contact)
        source_vol = float(vol[mask].sum())
        if source_vol <= 0:
            raise RuntimeError('Contact %d has no elements; refine the mesh' % contact)
        source = np.zeros(len(nodes))
        np.add.at(source, tets[mask].ravel(), np.repeat(vol[mask] / (4 * source_vol), 4))
        rhs = source - ret
        case = directory / ('contact%d' % contact)
        case.mkdir(exist_ok=True)
        results = case / 'results'
        text = sif(order, labels, sigma, contact, source_vol, area, pin, fibre, tags)
        if tree is None:
            tree = cKDTree(nodes * 1e-3)

        def check():
            v = read_potential(results, nodes, tree, resources.cpus)
            v -= v[pin]
            residual = action(nodes, tets, sigma, v, fibre=fibre) - rhs
            rel = float(np.linalg.norm(residual[free]) / max(np.linalg.norm(rhs[free]), 1e-300))
            return v, rel, float(residual[pin])

        checked = None
        if finished(case, text):
            # A job killed at its walltime loses only the contact it was on: an
            # earlier job's finished solve of this exact case is read back and
            # has to pass the same independent check as a fresh solve.
            resources.mark('verify contact %d (solved by an earlier job)' % contact)
            try:
                checked = check()
            except RuntimeError as error:
                print('Contact %d: earlier result unreadable (%s)' % (contact, error), flush=True)
            if checked and not accepted(*checked[1:]):
                print('Contact %d: earlier result fails the check (%.3e)' % (contact, checked[1]), flush=True)
                checked = None
        if checked is None:
            # Each solve gets fresh results; a failed rerun cannot import old VTUs.
            if results.exists():
                shutil.rmtree(results)
            results.mkdir()
            link = case / 'mesh'
            if not link.exists():
                link.symlink_to(directory.resolve() / 'mesh', target_is_directory=True)
            (case / 'case.sif').write_text(text)
            resources.mark('Elmer contact %d' % contact)
            command = [ELMER_BIN / 'ElmerSolver_mpi', 'case.sif']
            if resources.cpus > 1:
                command = [mpi, '-np', str(resources.cpus)] + command
            run(command, case, case / 'solver.log')
            resources.mark('verify contact %d' % contact)
            checked = check()
        v, rel, reaction = checked
        version = version or next((line.split('Version:', 1)[1].strip() for line
                                   in (case / 'solver.log').read_text(errors='replace').splitlines()
                                   if 'MAIN: Version:' in line), 'unknown')
        print('Contact %d: independently checked relative residual %.3e, pin reaction %.3e A' %
              (contact, rel, reaction), flush=True)
        if not accepted(rel, reaction):
            raise RuntimeError('Elmer contact %d failed independent residual/current checks; no solution published' % contact)
        phi[pos] = v
        residuals.append(rel); reactions.append(reaction)
        progress(round(12 + 83 * (pos + 1) / len(contact_ids)))
    V = np.zeros(len(nodes))
    if C.SOURCE_CONTACT in contact_ids and C.SINK_CONTACT in contact_ids:
        V = C.CURRENT_A * (phi[contact_ids.index(C.SOURCE_CONTACT)] - phi[contact_ids.index(C.SINK_CONTACT)])
    resources.mark('solution complete')
    extra = {} if fibre is None else dict(fibre_bin=fibre[0], fibre_a=fibre[1])
    return dict(extra, nodes=nodes, tets=tets, label=labels, sigma=sigma, phi=phi, V=V,
                relative_residuals=np.array(residuals), pin_reactions=np.array(reactions),
                accepted_residual=ACCEPTED_RESIDUAL, requested_residual=1e-11,
                contact_ids=np.array(contact_ids), order=np.array(order), pinned_node=pin,
                backend='Elmer StatCurrentSolver', elmer_bin=str(ELMER_BIN), elmer_version=version,
                cpus=resources.cpus,
                memory_budget_gb=resources.memory_gb, peak_rss_gb=resources.peak_gb)


def solve(nodes, tets, contact_stl, insulator_stl, include_background, directory, resources, progress):
    from assign_and_solve import classify_tets, sigma_of, dura_leak_report
    order = C.order_for(contact_stl)
    # Guard loading/classification, not the old all-at-once matrix estimate.
    minimum = (nodes.nbytes + tets.nbytes) / 1e9 * 3 + .5
    if minimum > resources.memory_gb:
        raise RuntimeError('Mesh preparation estimate %.2f GB exceeds %.2f GB local budget' % (minimum, resources.memory_gb))
    resources.mark('classification')
    lab = classify_tets(nodes, tets, contact_stl=contact_stl, insulator_stl=insulator_stl, order=order,
                        workers=resources.cpus)
    keep = np.ones(len(tets), bool) if include_background else lab >= 0
    active_tets, active_lab = tets[keep], lab[keep]
    joined, islands = drop_islands(nodes, active_tets, active_lab, order)
    if not joined.all():
        print('Dropped %d mesh pieces not connected to the lead: %d tets, %.2f mm3 (%s)'
              % (islands['components'] - islands['kept_components'], islands['dropped_tets'],
                 islands['dropped_mm3'], islands['dropped_by_tissue']), flush=True)
        active_tets, active_lab = active_tets[joined], active_lab[joined]
    del joined
    used = np.unique(active_tets)
    remap = np.full(len(nodes), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    active_tets = remap[active_tets]
    active_nodes = nodes[used]
    del keep, lab, remap, used
    print('Active: %d nodes, %d tetrahedra' % (len(active_nodes), len(active_tets)), flush=True)
    progress(5)
    nbr = face_neighbours(active_tets, resources.cpus)
    dura_leak_report(active_nodes, active_tets, active_lab, order=order, nbr=nbr)
    sigma = sigma_of(active_lab, include_background=include_background, order=order)
    resources.mark('fibre directions')
    fibre_bin, fibre_a, fibre_info = fibres.fibre_field(active_nodes, active_tets, active_lab, order)
    fibre = None if fibre_bin is None else (fibre_bin, fibre_a)
    print('White matter: %s' % (('anisotropic, %(bins_used)d of %(bins)d direction bins, '
                                 'cord turns %(turn_deg).1f deg, farthest element %(max_axis_distance_mm).1f mm '
                                 'from the axis' % fibre_info) if fibre else 'isotropic'), flush=True)
    result = solve_prepared(active_nodes, active_tets, active_lab, sigma, order,
                            sorted(contact_stl), directory, resources, progress, nbr=nbr, fibre=fibre)
    result['background_included'] = include_background
    import json
    result['fibre_info'] = json.dumps(fibre_info)
    result['dropped_islands'] = json.dumps(islands)
    return result
