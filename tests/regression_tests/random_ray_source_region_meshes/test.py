from glob import glob
import os
from pathlib import Path
import shutil

import h5py
import numpy as np
import openmc
import openmc.lib
from openmc.utility_funcs import change_directory
import pytest
from scipy.io import netcdf_file

from tests.regression_tests import config
from tests.testing_harness import TolerantPyAPITestHarness

UMESH_DIR = Path(__file__).parents[1] / 'unstructured_mesh'
FEATURES = {'xdg-moab': ['xdg', 'dagmc'], 'xdg-libmesh': ['xdg', 'libmesh'],
            'xdg-libmesh-unnormalized': ['xdg', 'libmesh']}
# Half of the width of the boxes that the models are a row of
W = 2.25


class SourceRegionMeshHarness(TolerantPyAPITestHarness):
    def _get_results(self):
        """The tallies, the weight windows, the volumes of unstructured meshes
        and the flux that libMesh writes"""
        lines = [super()._get_results()]
        lines += [f'{x:.6e}' for x in
                  openmc.WeightWindowsList.from_hdf5()[0].lower_ww_bounds.flat]
        with openmc.StatePoint(self._sp_name) as sp:
            lines += [f'{mesh.volumes.sum():.6e}'
                      for mesh in sp.meshes.values()
                      if isinstance(mesh, openmc.UnstructuredMesh)]
        for path in sorted(glob('tally_*.e')):
            # (NetCDF-4, which libMesh can be built to write, is HDF5)
            netcdf4 = h5py.is_hdf5(path)
            with (h5py.File(path, 'r') if netcdf4
                  else netcdf_file(path, mmap=False)) as exodus:
                variables = exodus if netcdf4 else exodus.variables
                names = [b''.join(name).decode().strip('\x00')
                         for name in variables['name_nod_var'][()]]
                for score in ('flux', 'total'):
                    index = names.index(f'{score}_total_mean') + 1
                    lines += [f'{x:.6e}' for x in
                              variables[f'vals_nod_var{index}'][-1]]
        return '\n'.join(lines)

    def _cleanup(self):
        super()._cleanup()
        for pattern in ('mgxs.h5', 'weight_windows.h5', '*tets.*',
                        'tally_*.e'):
            for path in glob(pattern):
                os.remove(path)


def write_tets(case):
    """Return the paths of the 24-tet mesh of the unit cube and of a copy with
    the node order of its tets flipped. For MOAB, the nodes in the middle of
    its faces move so that its volumes differ, and it is in two parts, so that
    rays leave it and come back."""
    if 'libmesh' in case:
        original = os.path.relpath(UMESH_DIR / 'test_mesh_dagmc_tets.e')
        shutil.copy(original, 'flipped_tets.e')
        with netcdf_file('flipped_tets.e', 'a', mmap=False) as exodus:
            connect = exodus.variables['connect1']
            connect[:] = connect[:][:, [1, 0, 2, 3]]
        return original, 'flipped_tets.e'

    lines = (UMESH_DIR / 'test_mesh_dagmc_tets.vtk').read_text().splitlines()
    i_points = next(i for i, line in enumerate(lines)
                    if line.startswith('POINTS'))
    n_points = int(lines[i_points].split()[1])
    points = np.loadtxt(lines[i_points + 1:i_points + 1 + n_points])
    points = np.where(points == 0.0, 0.4, points)
    points = np.vstack([points - [0.75, 0.0, 0.0], points + [0.75, 0.0, 0.0]])
    tets = np.array([line.split()[1:] for line in lines
                     if line.startswith('4 ')], dtype=int)
    tets = np.vstack([tets, tets + n_points])

    for path, nodes in (('tets.vtk', tets),
                        ('flipped_tets.vtk', tets[:, [1, 0, 2, 3]])):
        text = ['# vtk DataFile Version 2.0', 'tets', 'ASCII',
                'DATASET UNSTRUCTURED_GRID', f'POINTS {len(points)} double']
        text += [' '.join(map(str, point)) for point in points]
        text.append(f'CELLS {len(nodes)} {5 * len(nodes)}')
        text += ['4 ' + ' '.join(map(str, tet)) for tet in nodes]
        text.append(f'CELL_TYPES {len(nodes)}')
        text += ['10'] * len(nodes)
        Path(path).write_text('\n'.join(text) + '\n')
    return 'tets.vtk', 'flipped_tets.vtk'


def make_materials():
    """Write a one-group library and return its material and the materials"""
    groups = openmc.mgxs.EnergyGroups([1e-5, 20.0e6])
    xsdata = openmc.XSdata('mat', groups)
    xsdata.order = 0
    xsdata.set_total([0.2])
    xsdata.set_absorption([0.1])
    xsdata.set_scatter_matrix(np.array([[[0.1]]]))
    mgxs_lib = openmc.MGXSLibrary(groups)
    mgxs_lib.add_xsdata(xsdata)
    mgxs_lib.export_to_hdf5('mgxs.h5')
    mat = openmc.Material(name='mat')
    mat.set_density('macro', 1.0)
    mat.add_macroscopic(openmc.Macroscopic('mat'))
    materials = openmc.Materials([mat])
    materials.cross_sections = 'mgxs.h5'
    return mat, materials


def make_model(case):
    """A random ray model of a row of boxes along x, with a mesh over part of
    each, which rays leave and enter. The only source is a point outside of the
    meshes, with linear sources and weight windows on one of the meshes. A
    tally on a mesh counts only what that mesh subdivides, as
    otherwise the bin of a source region would depend on where a ray found it
    first."""
    mat, materials = make_materials()

    x_walls = ([-W, W, W + 1.0] if case.startswith('xdg')
               else [-W, W, 3.0 * W, 5.0 * W])
    # The last box is void, with a mesh nowhere near it
    x_walls.append(x_walls[-1] + 1.0)
    planes = [openmc.XPlane(x) for x in x_walls]
    planes[0].boundary_type = 'reflective'
    planes[-1].boundary_type = 'reflective'
    sides = -openmc.model.RectangularPrism(
        2.0 * W, 2.0 * W, axis='x', boundary_type='reflective')
    cells = [openmc.Cell(fill=mat, region=+a & -b & sides)
             for a, b in zip(planes[:-1], planes[1:])]
    cells[-1].fill = None
    far = openmc.RegularMesh()
    far.lower_left = (100.0, 100.0, 100.0)
    far.upper_right = (101.0, 101.0, 101.0)
    far.dimension = (2, 2, 2)

    if case == 'regular':
        # Part of a box, with volumes that differ
        part = openmc.RectilinearMesh()
        part.x_grid = [-0.5, 0.1, 0.5]
        part.y_grid = [-0.5, -0.2, 0.0, 0.5]
        part.z_grid = [-0.5, 0.4, 0.5]
        # All of a box but for a sliver on each side
        sliver = 2.5e-8
        most = openmc.RegularMesh()
        most.lower_left = (W + sliver, -W + sliver, -W + sliver)
        most.upper_right = (3.0 * W - sliver, W - sliver, W - sliver)
        most.dimension = (3, 3, 3)
        # Part of a box, in 2D
        flat = openmc.RegularMesh()
        flat.lower_left = (3.5 * W, -1.0)
        flat.upper_right = (4.5 * W, 1.0)
        flat.dimension = (2, 2)
        meshes = [(part, [cells[0]]), (most, [cells[1]]), (flat, [cells[2]])]
        mesh_filters = [[openmc.MeshFilter(mesh)] for mesh, _ in meshes]
        ww_mesh = part
    elif case == 'radial':
        # A hollow wedge, a cone, and part of a hollow sphere
        wedge = openmc.CylindricalMesh(
            r_grid=[0.4, 0.7, 1.0], z_grid=[-0.5, 0.0, 0.5],
            phi_grid=[1.0, 1.8, 2.5], origin=(0.3, -0.2, 0.1))
        cone = openmc.SphericalMesh(
            r_grid=[0.0, 0.5, 1.0], theta_grid=[0.0, 0.4, 0.8],
            phi_grid=[0.0, 2.0 * np.pi], origin=(2.0 * W, 0.2, -0.1))
        shell = openmc.SphericalMesh(
            r_grid=[0.4, 0.7, 1.0], theta_grid=[0.5, 1.5, 2.5],
            phi_grid=[0.0, np.pi, 2.0 * np.pi], origin=(4.0 * W, -0.1, 0.2))
        meshes = [(wedge, [cells[0]]), (cone, [cells[1]]), (shell, [cells[2]])]
        mesh_filters = [[openmc.MeshFilter(mesh)] for mesh, _ in meshes]
        ww_mesh = wedge
    else:
        path, flipped = write_tets(case)
        library = case.split('-')[1]

        def umesh(path, interface, multiplier=1.0):
            mesh = openmc.UnstructuredMesh(
                path, library=library, length_multiplier=multiplier)
            mesh.interface = interface
            return mesh

        # This mesh subdivides the first box. Tallies are on it through XDG,
        # and on it with its node order flipped through XDG and through the
        # native interface, which only looks up points, as do the weight
        # windows.
        main = umesh(path, 'xdg')
        ww_mesh = umesh(flipped, 'native')
        # This one only subdivides, and is scaled to cover the middle box. A
        # tally on that box has the same mesh with its node order flipped
        # through the native interface.
        scaled = umesh(path, 'xdg', multiplier=8.0)
        meshes = [(main, [cells[0]]), (scaled, [cells[1]])]
        mesh_filters = [
            [openmc.MeshFilter(main)],
            [openmc.MeshFilter(ww_mesh)],
            [openmc.MeshFilter(umesh(flipped, 'xdg'))],
            [openmc.CellFilter(cells[1]),
             openmc.MeshFilter(umesh(flipped, 'native', multiplier=8.0))],
        ]

    settings = openmc.Settings()
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    # (Linear sources start after the tenth batch)
    settings.batches = 15
    settings.inactive = 5
    settings.particles = 100
    settings.random_ray['distance_active'] = 30.0
    settings.random_ray['distance_inactive'] = 3.0
    settings.random_ray['ray_source'] = openmc.IndependentSource(
        space=openmc.stats.Box((x_walls[0], -W, -W), (x_walls[-1], W, W)))
    settings.random_ray['source_shape'] = 'linear'
    # (Weight windows and Exodus output take either kind of flux tally)
    settings.random_ray['volume_normalized_flux_tallies'] = (
        case not in ('xdg-moab', 'xdg-libmesh-unnormalized'))
    settings.random_ray['source_region_meshes'] = meshes + [(far, [cells[-1]])]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((-2.0, 2.0, 2.0)),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))
    generator = openmc.WeightWindowGenerator(
        mesh=ww_mesh, max_realizations=5,
        method='fw_cadis' if case == 'xdg-libmesh' else 'magic')
    if case == 'radial':
        generator.update_parameters = {
            'value': 'rel_err', 'threshold': 1.0, 'ratio': 5.0}
    settings.weight_window_generators = generator

    tallies = openmc.Tallies()
    for filters in [[openmc.CellFilter(cells)]] + mesh_filters:
        tally = openmc.Tally()
        tally.filters = filters
        tally.scores = ['flux', 'total']
        tally.estimator = 'tracklength'
        tallies.append(tally)

    return openmc.Model(openmc.Geometry(cells), materials, settings, tallies)


@pytest.mark.parametrize(
    "case", ['regular', 'radial', 'xdg-moab', 'xdg-libmesh',
             'xdg-libmesh-unnormalized'])
def test_random_ray_source_region_meshes(case, capfd):
    """Source region meshes that cover only part of their cells, with tallies
    and weight windows on them. What each part of the model guards:

    - Meshes over part of a box, in two parts, or short of it by a sliver:
      what is outside of a mesh is attenuated, as one more bin, in its place
      along the ray.
    - The point source outside of the meshes: it is not dropped.
    - The void box with a mesh nowhere near it: its flux is tallied, and the
      mesh is warned about, once.
    - radial: cylindrical and spherical meshes that do not span all angles
      find tracks that start outside of them.
    - xdg-*: subdivision through XDG, with mesh tallies through XDG and
      through the native interface, and weight windows through the native
      one.
    - Volumes that differ and flipped node order: weight windows and Exodus
      output are divided by volume once, and volumes are positive.
    - The cases differ in whether flux tallies are normalized by volume, and
      xdg-libmesh generates FW-CADIS weight windows, which writes the forward
      Exodus file.
    """
    if not all(openmc.lib.feature_enabled(f) for f in FEATURES.get(case, [])):
        pytest.skip(f"{case} is not enabled.")
    # A ray handed to another rank restarts on a bin boundary. On a radial
    # surface of these meshes, going inward at a grazing angle, it misses its
    # way out again, as find_r_crossing takes it to still be on the surface.
    if config['mpi'] and case == 'radial':
        pytest.skip("Radial crossings of restarted rays")

    with change_directory(case):
        openmc.reset_auto_ids()
        model = make_model(case)
        harness = SourceRegionMeshHarness('statepoint.15.h5', model)
        harness.main()

    # Only the mesh nowhere near its cell is reported, once
    if not config['build_inputs']:
        far = model.settings.random_ray['source_region_meshes'][-1][0]
        output = ' '.join(capfd.readouterr().out.split())
        assert output.count('so it subdivides nothing') == 1
        assert f'source region mesh {far.id} is applied' in output


@pytest.mark.parametrize("library,use", [
    ('moab', 'subdivision'),
    ('libmesh', 'subdivision'),
    ('libmesh', 'tally'),
])
def test_native_interface_errors(library, use, tmp_path):
    """What the native interfaces leave to XDG is reported up front"""
    if not openmc.lib.feature_enabled(
            'dagmc' if library == 'moab' else 'libmesh'):
        pytest.skip(f"{library} is not enabled.")
    if config['build_inputs']:
        return

    openmc.reset_auto_ids()
    extension = 'vtk' if library == 'moab' else 'e'
    mesh = openmc.UnstructuredMesh(
        UMESH_DIR / f'test_mesh_dagmc_tets.{extension}', library=library)
    with change_directory(tmp_path):
        model = make_model('regular')
    model.settings.weight_window_generators = []
    if use == 'subdivision':
        model.settings.random_ray['source_region_meshes'] = [
            (mesh, [model.geometry.root_universe])]
    else:
        model.settings.random_ray = {}
        model.tallies[0].filters = [openmc.MeshFilter(mesh)]

    # (Aborting can lose the message with many ranks)
    if config['mpi']:
        mpi_args, message = [config['mpiexec'], '-n', config['mpi_np']], None
    else:
        mpi_args, message = None, 'XDG interface'
    with pytest.raises(RuntimeError, match=message):
        model.run(cwd=tmp_path, openmc_exec=config['exe'], mpi_args=mpi_args)


def make_hand_off_model():
    """A lattice whose plane lies just before the surface of its cell, with a
    mesh plane between the two. A ray through the first tile heads for that
    surface, as the lattice plane is closer by less than the tolerance of the
    distance to a boundary, and so runs past the lattice plane, where the mesh
    trace gives it a piece. With two ranks, the subdomains meet at x = 0 until
    they are rebalanced after the first batch, so that piece, at x = 0.5, is on
    the other rank's side."""
    mat, materials = make_materials()
    gap = 1.0e-6
    tile = openmc.Universe(cells=[openmc.Cell(fill=mat)])
    lattice = openmc.RectLattice()
    lattice.lower_left = (-2.0, -0.5, -0.5)
    lattice.pitch = (2.5 - gap, 1.0, 1.0)
    lattice.universes = [[[tile, tile]]]
    planes = [openmc.XPlane(-2.0, boundary_type='reflective'),
              openmc.XPlane(0.5),
              openmc.XPlane(2.0, boundary_type='reflective')]
    sides = -openmc.model.RectangularPrism(
        1.0, 1.0, axis='x', boundary_type='reflective')
    cells = [openmc.Cell(fill=lattice, region=+planes[0] & -planes[1] & sides),
             openmc.Cell(fill=mat, region=+planes[1] & -planes[2] & sides)]
    geometry = openmc.Geometry(cells)

    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [-2.0, 0.5 - gap / 2.0, 2.0]
    mesh.y_grid = [-0.5, 0.5]
    mesh.z_grid = [-0.5, 0.5]

    settings = openmc.Settings()
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    settings.batches = 1
    settings.inactive = 0
    settings.particles = 30
    # (One lost ray stops the run)
    settings.max_lost_particles = 1
    settings.random_ray['distance_active'] = 20.0
    settings.random_ray['distance_inactive'] = 4.0
    settings.random_ray['ray_source'] = openmc.IndependentSource(
        space=openmc.stats.Box((-2.0, -0.5, -0.5), (2.0, 0.5, 0.5)))
    settings.random_ray['source_region_meshes'] = [
        (mesh, [geometry.root_universe])]
    settings.source = openmc.IndependentSource(
        energy=openmc.stats.Discrete([1.0e6], [1.0]),
        constraints={'domains': [mat]})

    return openmc.Model(geometry, materials, settings)


def test_hand_off_past_lattice_plane(tmp_path):
    """A ray is not handed to another rank for a piece of its segment past a
    lattice plane, from where it would be lost"""
    if not config['mpi']:
        pytest.skip("Domain decomposition needs MPI.")
    if config['build_inputs']:
        return

    openmc.reset_auto_ids()
    with change_directory(tmp_path):
        model = make_hand_off_model()
    model.run(cwd=tmp_path, openmc_exec=config['exe'],
              mpi_args=[config['mpiexec'], '-n', config['mpi_np']])
