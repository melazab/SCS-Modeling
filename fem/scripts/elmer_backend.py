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
from element_batches import elements, action, batch_size
from mesh_faces import face_batches

ELMER_BIN = Path('/opt/elmerfem/bin')


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
    """Reject disconnected active domains using a batched boolean graph."""
    graph = sp.csr_matrix((n, n), dtype=bool)
    for start in range(0, len(tets), batch_size()):
        tet = tets[start:start + batch_size()]
        row = np.repeat(tet[:, 0], 3)
        col = tet[:, 1:].ravel()
        graph += sp.coo_matrix((np.ones(len(row), bool), (row, col)), shape=(n, n)).tocsr()
    count = connected_components(graph, directed=False, return_labels=False)
    if count != 1:
        raise RuntimeError('Active mesh has %d disconnected components; check tissue connectivity before solving.' % count)


def write_mesh(directory, nodes, tets, labels, order, pin=None):
    """Write native Elmer mesh in SI units; return normalized nodal return."""
    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / 'mesh.nodes', 'w') as stream:
        for start in range(0, len(nodes), batch_size()):
            xyz = nodes[start:start + batch_size()] * 1e-3
            ids = np.arange(start + 1, start + len(xyz) + 1)
            np.savetxt(stream, np.column_stack((ids, np.full(len(ids), -1), xyz)),
                       fmt=['%d', '%d', '%.17g', '%.17g', '%.17g'])
    with open(directory / 'mesh.elements', 'w') as stream:
        for start in range(0, len(tets), batch_size()):
            tet = tets[start:start + batch_size()]
            lab = labels[start:start + len(tet)]
            tag = np.where(lab >= 0, lab + 1, len(order) + 1)
            ids = np.arange(start + 1, start + len(tet) + 1)
            np.savetxt(stream, np.column_stack((ids, tag, np.full(len(tet), 504), tet + 1)), fmt='%d')
    nodal_area = np.zeros(len(nodes))
    if pin is None:
        pin = int(np.argmax(np.linalg.norm(nodes[:, :2] - np.array([55.60, 87.0]), axis=1)))
    count = 0
    with open(directory / 'mesh.boundary', 'w') as stream:
        for faces, owner, paired, lone in face_batches(tets):
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


def sif(order, labels, sigma, contact, source_volume, return_area, pin):
    bodies = []
    for label in np.unique(labels):
        tag = int(label) + 1 if label >= 0 else len(order) + 1
        name = order[int(label)] if label >= 0 else 'background'
        cond = float(sigma[np.flatnonzero(labels == label)[0]])
        bodyforce = '\n  Body Force = 1' if name == 'contact%d' % contact else ''
        bodies.append('Body %d\n  Target Bodies(1) = %d\n  Equation = 1\n  Material = %d%s\nEnd\n'
                      'Material %d\n  Electric Conductivity = %.17g\nEnd\n' %
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


def read_potential(results, nodes):
    """Map partition output to original nodes and check interface agreement."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    tree = cKDTree(nodes * 1e-3)
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
        dist, ids = tree.query(vtk_to_numpy(grid.GetPoints().GetData()))
        if not np.isfinite(v).all() or np.max(dist) > 1e-9:
            raise RuntimeError('Elmer result coordinates/values do not match input mesh')
        present = np.isfinite(values[ids])
        if present.any() and not np.allclose(values[ids[present]], v[present], rtol=1e-8, atol=1e-8):
            raise RuntimeError('MPI partitions disagree at shared nodes')
        values[ids] = v
    if not np.isfinite(values).all():
        raise RuntimeError('Elmer result is missing nodes')
    return values


def solve_prepared(nodes, tets, labels, sigma, order, contact_ids, directory, resources, progress=lambda p: None):
    """Solve a classified mesh; independently verify every unit-current field."""
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
    ret, area = write_mesh(directory / 'mesh', nodes, tets, labels, order, pin=pin)
    progress(12)
    if resources.cpus > 1:
        resources.mark('mesh partitioning')
        run([ELMER_BIN / 'ElmerGrid', '2', '2', 'mesh', '-metis', str(resources.cpus), '-parthypre'],
            directory, directory / 'partition.log')
    phi = np.empty((len(contact_ids), len(nodes)))
    residuals = []
    reactions = []
    free = np.ones(len(nodes), bool); free[pin] = False
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
        # Each solve gets fresh results; a failed rerun cannot import old VTUs.
        results = case / 'results'
        if results.exists():
            shutil.rmtree(results)
        results.mkdir()
        link = case / 'mesh'
        if not link.exists():
            link.symlink_to(directory.resolve() / 'mesh', target_is_directory=True)
        (case / 'case.sif').write_text(sif(order, labels, sigma, contact, source_vol, area, pin))
        resources.mark('Elmer contact %d' % contact)
        command = [ELMER_BIN / 'ElmerSolver_mpi', 'case.sif']
        if resources.cpus > 1:
            command = [mpi, '-np', str(resources.cpus)] + command
        run(command, case, case / 'solver.log')
        resources.mark('verify contact %d' % contact)
        v = read_potential(results, nodes)
        v -= v[pin]
        residual = action(nodes, tets, sigma, v) - rhs
        rel = float(np.linalg.norm(residual[free]) / max(np.linalg.norm(rhs[free]), 1e-300))
        reaction = float(residual[pin])
        print('Contact %d: independently checked relative residual %.3e, pin reaction %.3e A' %
              (contact, rel, reaction), flush=True)
        if not np.isfinite(rel) or rel > ACCEPTED_RESIDUAL or abs(reaction) > 1e-6:
            raise RuntimeError('Elmer contact %d failed independent residual/current checks; no solution published' % contact)
        phi[pos] = v
        residuals.append(rel); reactions.append(reaction)
        progress(round(12 + 83 * (pos + 1) / len(contact_ids)))
    V = np.zeros(len(nodes))
    if C.SOURCE_CONTACT in contact_ids and C.SINK_CONTACT in contact_ids:
        V = C.CURRENT_A * (phi[contact_ids.index(C.SOURCE_CONTACT)] - phi[contact_ids.index(C.SINK_CONTACT)])
    resources.mark('solution complete')
    return dict(nodes=nodes, tets=tets, label=labels, sigma=sigma, phi=phi, V=V,
                relative_residuals=np.array(residuals), pin_reactions=np.array(reactions),
                accepted_residual=ACCEPTED_RESIDUAL, requested_residual=1e-11,
                contact_ids=np.array(contact_ids), order=np.array(order), pinned_node=pin,
                backend='Elmer StatCurrentSolver', cpus=resources.cpus,
                memory_budget_gb=resources.memory_gb, peak_rss_gb=resources.peak_gb)


def solve(nodes, tets, contact_stl, insulator_stl, include_background, directory, resources, progress):
    from assign_and_solve import classify_tets, sigma_of, dura_leak_report
    order = C.order_for(contact_stl)
    # Guard loading/classification, not the old all-at-once matrix estimate.
    minimum = (nodes.nbytes + tets.nbytes) / 1e9 * 3 + .5
    if minimum > resources.memory_gb:
        raise RuntimeError('Mesh preparation estimate %.2f GB exceeds %.2f GB local budget' % (minimum, resources.memory_gb))
    resources.mark('classification')
    lab = classify_tets(nodes, tets, contact_stl=contact_stl, insulator_stl=insulator_stl, order=order)
    keep = np.ones(len(tets), bool) if include_background else lab >= 0
    active_tets, active_lab = tets[keep], lab[keep]
    used = np.unique(active_tets)
    remap = np.full(len(nodes), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    active_tets = remap[active_tets]
    active_nodes = nodes[used]
    del keep, lab, remap, used
    print('Active: %d nodes, %d tetrahedra' % (len(active_nodes), len(active_tets)), flush=True)
    progress(5)
    dura_leak_report(active_nodes, active_tets, active_lab, order=order)
    sigma = sigma_of(active_lab, include_background=include_background, order=order)
    result = solve_prepared(active_nodes, active_tets, active_lab, sigma, order,
                            sorted(contact_stl), directory, resources, progress)
    result['background_included'] = include_background
    return result
