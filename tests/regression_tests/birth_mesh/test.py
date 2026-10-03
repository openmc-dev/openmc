"""Test of the fission bank birth-mesh tags in a decoupled two-slab problem.

Two one-group fissile slabs are separated by a 10 cm pure absorber with
Sigma_t = 10/cm, in which a neutron is absorbed at its first collision. A
neutron needs at least 10 cm of uncollided path in the absorber to cross it,
so no neutron crosses it except with probability <= e^-100 (and 2^-64 per
sampled flight distance, for a random number of exactly 0, which gives an
infinite flight). The birth position of the neutron that produced a banked
fission site is therefore in the same slab as the site itself. With a birth
mesh that has one bin per slab, every tag must equal the bin of the position
of its own site, so the test does not need reference results. A finer mesh,
with two bins per slab, checks that the tags are taken from the birth position
of the neutron that produced the site rather than from the position of the
site, for a regular and for a nonuniform rectilinear mesh.

The tests also check which output files carry the tags, the restart path that
reads them back, and the rejection of a tag dataset whose length differs from
that of the source bank.

"""
import shutil

import h5py
import numpy as np
import openmc
import openmc.lib
import pytest

from tests.regression_tests import config

# One-group data: subcritical fuel, so that the fixed-source variant with
# fission neutrons is well defined, and a pure absorber
FUEL = {'total': 1.0, 'absorption': 0.5, 'fission': 0.15, 'nu': 2.5}
ABSORBER_SIGMA = 10.0  # 1/cm, total = absorption
HALF_WIDTH = 5.0  # cm, half-thickness of the absorber
FUEL_WIDTH = 5.0  # cm, thickness of each fuel slab


def make_model(library):
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    fuel_xs = openmc.XSdata('fuel', groups)
    fuel_xs.order = 0
    fuel_xs.set_total([FUEL['total']])
    fuel_xs.set_absorption([FUEL['absorption']])
    fuel_xs.set_fission([FUEL['fission']])
    fuel_xs.set_nu_fission([FUEL['nu'] * FUEL['fission']])
    fuel_xs.set_scatter_matrix([[[FUEL['total'] - FUEL['absorption']]]])
    fuel_xs.set_chi([1.0])
    absorber_xs = openmc.XSdata('absorber', groups)
    absorber_xs.order = 0
    absorber_xs.set_total([ABSORBER_SIGMA])
    absorber_xs.set_absorption([ABSORBER_SIGMA])
    absorber_xs.set_scatter_matrix([[[0.0]]])
    mgxs_library = openmc.MGXSLibrary(groups)
    mgxs_library.add_xsdata(fuel_xs)
    mgxs_library.add_xsdata(absorber_xs)
    mgxs_library.export_to_hdf5(library)

    fuel = openmc.Material(name='fuel')
    fuel.set_density('macro', 1.0)
    fuel.add_macroscopic('fuel')
    absorber = openmc.Material(name='absorber')
    absorber.set_density('macro', 1.0)
    absorber.add_macroscopic('absorber')

    x_max = HALF_WIDTH + FUEL_WIDTH
    x0 = openmc.XPlane(-x_max, boundary_type='reflective')
    x1 = openmc.XPlane(-HALF_WIDTH)
    x2 = openmc.XPlane(HALF_WIDTH)
    x3 = openmc.XPlane(x_max, boundary_type='reflective')
    y0 = openmc.YPlane(-5.0, boundary_type='reflective')
    y1 = openmc.YPlane(5.0, boundary_type='reflective')
    z0 = openmc.ZPlane(-5.0, boundary_type='reflective')
    z1 = openmc.ZPlane(5.0, boundary_type='reflective')
    yz = +y0 & -y1 & +z0 & -z1
    cells = [openmc.Cell(fill=fuel, region=+x0 & -x1 & yz),
             openmc.Cell(fill=absorber, region=+x1 & -x2 & yz),
             openmc.Cell(fill=fuel, region=+x2 & -x3 & yz)]

    model = openmc.Model()
    model.materials = openmc.Materials([fuel, absorber])
    model.materials.cross_sections = str(library)
    model.geometry = openmc.Geometry(cells)
    model.settings.energy_mode = 'multi-group'
    model.settings.particles = 2000
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-x_max, -5.0, -5.0), (x_max, 5.0, 5.0)),
        constraints={'fissionable': True})
    return model, x1


def slab_mesh(nx):
    """Regular mesh over the problem with nx bins of equal width along x."""
    x_max = HALF_WIDTH + FUEL_WIDTH
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-x_max, -5.0, -5.0)
    mesh.upper_right = (x_max, 5.0, 5.0)
    mesh.dimension = (nx, 1, 1)
    return mesh


def rectilinear_slab_mesh():
    """Nonuniform rectilinear mesh with two bins along x in each fuel slab, one
    for the absorber, and two unequal bins along y."""
    x_max = HALF_WIDTH + FUEL_WIDTH
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [-x_max, -7.0, -HALF_WIDTH, HALF_WIDTH, 7.5, x_max]
    mesh.y_grid = [-5.0, 1.0, 5.0]
    mesh.z_grid = [-5.0, 5.0]
    return mesh


def bin_of(mesh, r):
    """0-based bin (x fastest) of each position of r in a Cartesian mesh."""
    xyz = np.column_stack([r['x'], r['y'], r['z']])
    if isinstance(mesh, openmc.RegularMesh):
        ll = np.asarray(mesh.lower_left)
        width = (np.asarray(mesh.upper_right) - ll) / mesh.dimension
        ijk = np.floor((xyz - ll) / width).astype(int)
    else:
        grids = (mesh.x_grid, mesh.y_grid, mesh.z_grid)
        ijk = np.column_stack([np.searchsorted(g, xyz[:, i], side='right') - 1
                               for i, g in enumerate(grids)])
    nx, ny, _ = mesh.dimension
    return ijk[:, 0] + nx * (ijk[:, 1] + ny * ijk[:, 2])


def run_model(model, path, **kwargs):
    kwargs.update({'cwd': path, 'openmc_exec': config['exe'],
                   'event_based': config['event']})
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    return model.run(**kwargs)


def eigenvalue_model(tmp_path, mesh):
    model, surface = make_model(tmp_path / 'two_slab.h5')
    model.settings.run_mode = 'eigenvalue'
    model.settings.batches = 6
    model.settings.inactive = 2
    model.settings.birth_mesh = mesh
    return model, surface


def read_source(path):
    with h5py.File(path, 'r') as f:
        attrs = dict(f.attrs)
        source = f['source_bank'][()]
        tags = f['birth_mesh_bin'][()] if 'birth_mesh_bin' in f else None
    return attrs, source, tags


def test_birth_mesh_eigenvalue(tmp_path):
    # One bin per slab: bin 0 holds the left slab and bin 1 the right one
    mesh = slab_mesh(2)
    model, surface = eigenvalue_model(tmp_path, mesh)
    # Also write a surface source file and the initial source file
    model.settings.surf_source_write = {'surface_ids': [surface.id],
                                        'max_particles': 1000}
    model.settings.write_initial_source = True
    sp_path = run_model(model, tmp_path)

    attrs, source, tags = read_source(sp_path)
    with openmc.StatePoint(sp_path) as sp:
        assert sp.source_present
        assert sp.birth_mesh_id == mesh.id
        np.testing.assert_array_equal(sp.birth_mesh_bin, tags)

    # The tags are int32, aligned with the source bank, and refer to the mesh
    assert attrs['birth_mesh_id'] == mesh.id
    assert tags is not None
    assert tags.dtype == np.int32
    assert tags.shape == source.shape == (model.settings.particles,)

    # Every site is tagged with the bin of its own slab, and both slabs hold
    # sites
    np.testing.assert_array_equal(tags, bin_of(mesh, source['r']))
    assert set(np.unique(tags)) == {0, 1}

    # Surface source and initial source sites are not fission bank sites, so
    # these files carry no tags
    for name in ('surface_source.h5', 'initial_source.h5'):
        file_attrs, _, file_tags = read_source(tmp_path / name)
        assert 'birth_mesh_id' not in file_attrs
        assert file_tags is None


def test_birth_mesh_source_files(tmp_path):
    mesh = slab_mesh(2)
    model, _ = eigenvalue_model(tmp_path, mesh)
    # Write the source bank of the last batch to source.6.h5 and the latest
    # source bank to source.h5 instead of to the state point
    model.settings.sourcepoint = {'separate': True, 'overwrite': True}
    sp_path = run_model(model, tmp_path)

    # The state point records the mesh but holds no source bank or tags
    with h5py.File(sp_path, 'r') as f:
        assert f.attrs['birth_mesh_id'] == mesh.id
        assert 'source_bank' not in f
        assert 'birth_mesh_bin' not in f
    with openmc.StatePoint(sp_path) as sp:
        assert not sp.source_present
        assert sp.birth_mesh_id == mesh.id
        assert sp.birth_mesh_bin is None

    # Both source point files carry the tags of the final source bank
    _, source, tags = read_source(tmp_path / 'source.6.h5')
    np.testing.assert_array_equal(tags, bin_of(mesh, source['r']))
    for name in ('source.6.h5', 'source.h5'):
        file_attrs, file_source, file_tags = read_source(tmp_path / name)
        assert file_attrs['birth_mesh_id'] == mesh.id
        assert file_tags.dtype == np.int32
        np.testing.assert_array_equal(file_source, source)
        np.testing.assert_array_equal(file_tags, tags)


@pytest.mark.parametrize('rectilinear', [False, True])
def test_birth_mesh_birth_position(tmp_path, rectilinear):
    # With two bins along x in each fuel slab (the first two and the last two),
    # a neutron is often born in a different bin from the one in which it
    # causes a fission, but always in the same slab
    mesh = rectilinear_slab_mesh() if rectilinear else slab_mesh(8)
    model, _ = eigenvalue_model(tmp_path, mesh)
    _, source, tags = read_source(run_model(model, tmp_path))

    nx = mesh.dimension[0]
    own_bin = bin_of(mesh, source['r'])
    assert tags.min() >= 0
    left, right = {0, 1}, {nx - 2, nx - 1}
    assert set(np.unique(own_bin % nx)) == left | right
    for slab in (left, right):
        in_slab = np.isin(own_bin % nx, list(slab))
        assert set(np.unique(tags[in_slab] % nx)) <= slab
    assert set(np.unique(tags // nx)) == set(np.unique(own_bin // nx))
    assert np.count_nonzero(tags != own_bin) > 0


def test_birth_mesh_restart(tmp_path, monkeypatch):
    mesh = slab_mesh(2)
    model, _ = eigenvalue_model(tmp_path, mesh)
    model.settings.statepoint = {'batches': [3, 6]}
    run_model(model, tmp_path)
    restart_file = tmp_path / 'statepoint.3.h5'
    _, source, tags = read_source(restart_file)

    # A restart run reads the tagged source bank back and tags the new sites
    sp_path = run_model(model, tmp_path, restart_file=restart_file)
    _, new_source, new_tags = read_source(sp_path)
    np.testing.assert_array_equal(new_tags, bin_of(mesh, new_source['r']))

    # The tags read from the restart file are those of the source sites
    if not config['mpi']:
        monkeypatch.chdir(tmp_path)
        model.init_lib(restart_file=restart_file, output=False)
        try:
            openmc.lib.simulation_init()
            bank = openmc.lib.source_bank()
            np.testing.assert_array_equal(bank['birth_mesh_bin'], tags)
            np.testing.assert_array_equal(bank['r'][:, 0], source['r']['x'])
            openmc.lib.simulation_finalize()
        finally:
            model.finalize_lib()

    # A tag dataset whose length differs from that of the source bank is
    # rejected, both on restart and for a file source. Without the length
    # check, a dataset that is one entry too long would be read silently on
    # restart, where each process reads only its share of the sites.
    length_error = ('birth_mesh_bin dataset in the source file must be '
                    'one-dimensional with the same length as the source_bank')
    initial_source = model.settings.source
    for n_tags in (tags.size - 1, tags.size + 1):
        bad_restart = tmp_path / f'statepoint.bad{n_tags}.h5'
        shutil.copy(restart_file, bad_restart)
        with h5py.File(bad_restart, 'r+') as f:
            del f['birth_mesh_bin']
            f.create_dataset('birth_mesh_bin',
                             data=np.zeros(n_tags, dtype=np.int32))
        with pytest.raises(RuntimeError, match=length_error):
            run_model(model, tmp_path, restart_file=bad_restart)

        bad_source = tmp_path / f'source.bad{n_tags}.h5'
        with h5py.File(bad_source, 'w') as f:
            f.attrs['filetype'] = np.bytes_('source')
            with h5py.File(restart_file, 'r') as sp:
                f.attrs['version'] = sp.attrs['version']
            f.create_dataset('source_bank', data=source)
            f.create_dataset('birth_mesh_bin',
                             data=np.zeros(n_tags, dtype=np.int32))
        model.settings.source = openmc.FileSource(bad_source)
        with pytest.raises(RuntimeError, match=length_error):
            run_model(model, tmp_path)

        # Restore the initial source so that the next restart check reads
        # only the restart file
        model.settings.source = initial_source


def test_birth_mesh_fixed_source(tmp_path, capsys):
    mesh = slab_mesh(2)
    model, surface = make_model(tmp_path / 'two_slab.h5')
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 2
    model.settings.birth_mesh = mesh
    model.settings.surf_source_write = {'surface_ids': [surface.id],
                                        'max_particles': 1000}
    sp_path = run_model(model, tmp_path)

    # The element is ignored with a warning, and no file carries tags
    assert 'no effect in fixed-source mode' in capsys.readouterr().out
    with h5py.File(sp_path, 'r') as f:
        assert 'birth_mesh_id' not in f.attrs
        assert 'birth_mesh_bin' not in f
    with openmc.StatePoint(sp_path) as sp:
        assert sp.birth_mesh_id is None
        assert sp.birth_mesh_bin is None
    file_attrs, _, file_tags = read_source(tmp_path / 'surface_source.h5')
    assert 'birth_mesh_id' not in file_attrs
    assert file_tags is None
