"""Tests of source files written and read by multigroup simulations.

In a multigroup simulation the source bank stores the energy group index of
each site in E. The source point files written from it are marked with an
energy_mode attribute equal to "multi-group", and a FileSource of such a file,
or of a multigroup state point, keeps the stored group of each site instead of
converting E from eV. Files without the marker, such as surface source files
and files written by openmc.write_source_file, hold energies in eV.

The tests that run OpenMC in a subprocess come first; the tests that sample
the source through openmc.lib follow.
"""
import math
from pathlib import Path
import re
import shutil
import sys

import h5py
import numpy as np
import openmc
import openmc.lib
import pytest

from tests.regression_tests import config

# Group edges in [eV] and fission spectrum of each library. Group 0 is the
# highest-energy group. The medium is fissile only in the group that the
# fission spectrum feeds, so every fission site is created in that group.
LIBRARIES = {
    '1g': ([1.0e-5, 2.0e7], [1.0]),
    '1g-zero-edge': ([0.0, 2.0e7], [1.0]),
    '2g-thermal': ([1.0e-5, 0.625, 2.0e7], [0.0, 1.0]),
    '2g-fast': ([1.0e-5, 0.625, 2.0e7], [1.0, 0.0]),
    '2g-fast-zero-edge': ([0.0, 0.625, 2.0e7], [1.0, 0.0]),
}

# Energy bounds wider than every group structure above
WIDE_BOUNDS = (-1.0, 3.0e7)

# Average energies of the two groups of the two-group libraries
LOW_AVERAGE = 0.5*(1.0e-5 + 0.625)
HIGH_AVERAGE = 0.5*(0.625 + 2.0e7)


def birth_group(name):
    """Energy group in which the fission neutrons of a library are born"""
    return LIBRARIES[name][1].index(1.0)


def make_model(name, directory='.', plane=None):
    """Reflective 10 cm cube of one material, with its library

    The cube is divided into an inner cell [-2, 2]^3 and an outer cell or, if
    a plane is given, into the two half-spaces of that plane. Each group has a
    total cross section of 1 /cm, half of it absorption and half scattering
    within the group, so particles never change group. The library is written
    to mgxs.h5 in the given directory. Returns the model and its two cells.

    """
    edges, chi = LIBRARIES[name]
    groups = openmc.mgxs.EnergyGroups(edges)
    n = groups.num_groups
    fission = 0.2*np.array(chi)
    xsdata = openmc.XSdata('medium', groups)
    xsdata.order = 0
    xsdata.set_total(np.ones(n))
    xsdata.set_absorption(np.full(n, 0.5))
    xsdata.set_scatter_matrix(np.diag(np.full(n, 0.5))[..., np.newaxis])
    xsdata.set_fission(fission)
    xsdata.set_nu_fission(2.0*fission)
    xsdata.set_chi(chi)
    library = openmc.MGXSLibrary(groups)
    library.add_xsdata(xsdata)
    library.export_to_hdf5(Path(directory) / 'mgxs.h5')

    mat = openmc.Material()
    mat.set_density('macro', 1.0)
    mat.add_macroscopic('medium')
    cube = openmc.model.RectangularParallelepiped(
        -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, boundary_type='reflective')
    if plane is None:
        inner = openmc.model.RectangularParallelepiped(
            -2.0, 2.0, -2.0, 2.0, -2.0, 2.0)
        cells = [openmc.Cell(fill=mat, region=-inner),
                 openmc.Cell(fill=mat, region=+inner & -cube)]
    else:
        cells = [openmc.Cell(fill=mat, region=-plane & -cube),
                 openmc.Cell(fill=mat, region=+plane & -cube)]

    model = openmc.Model()
    model.materials = openmc.Materials([mat])
    model.materials.cross_sections = 'mgxs.h5'
    model.geometry = openmc.Geometry(cells)
    model.settings.energy_mode = 'multi-group'
    model.settings.particles = 1000
    model.settings.batches = 3
    # Histories here take tens of events; the limit stops one that would run
    # forever, such as a particle in an invalid group if a test fails
    model.settings.max_particle_events = 1000
    return model, cells


def box_source(energies, probabilities):
    return openmc.IndependentSource(
        space=openmc.stats.Box((-5.0, -5.0, -5.0), (5.0, 5.0, 5.0)),
        energy=openmc.stats.Discrete(energies, probabilities))


def add_flux_tally(model, name):
    tally = openmc.Tally()
    tally.filters = [openmc.EnergyFilter(LIBRARIES[name][0])]
    tally.scores = ['flux']
    model.tallies = [tally]
    return tally


def run(model, cwd='.', mpi=False, **kwargs):
    """Run OpenMC in a subprocess, with MPI if requested and enabled"""
    if mpi and config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    kwargs.setdefault('output', False)
    return model.run(cwd=cwd, openmc_exec=config['exe'],
                     event_based=config['event'], **kwargs)


def run_flux(model, name, **kwargs):
    """Run the model and return its flux in each energy bin, lowest first"""
    tally = add_flux_tally(model, name)
    sp_path = run(model, **kwargs)
    with openmc.StatePoint(sp_path) as sp:
        return sp.tallies[tally.id].mean.ravel()


def write_source_points(directory, name, statepoint=False, mpi=False):
    """Run a multigroup eigenvalue simulation that writes source points

    By default the sites of the last batch are written to source.3.h5 and
    source.h5; with statepoint=True they are kept in statepoint.3.h5 instead.

    """
    model, _ = make_model(name, directory)
    model.settings.run_mode = 'eigenvalue'
    model.settings.inactive = 1
    # Start the particles in the group in which the medium is fissile
    model.settings.source = box_source(
        [0.1 if birth_group(name) == 1 else 1.0e6], [1.0])
    if statepoint:
        model.settings.sourcepoint = {'batches': [3]}
    else:
        model.settings.sourcepoint = {
            'batches': [3], 'separate': True, 'overwrite': True}
    run(model, cwd=directory, mpi=mpi)


@pytest.fixture(scope='module')
def mg_files(tmp_path_factory):
    """Directory with the source points of a library, written once"""
    directories = {}

    def get(name, statepoint=False):
        if (name, statepoint) not in directories:
            directory = tmp_path_factory.mktemp('mg_source')
            write_source_points(directory, name, statepoint)
            directories[name, statepoint] = directory
        return directories[name, statepoint]

    return get


@pytest.fixture(scope='module')
def surface_source_file(tmp_path_factory):
    """Surface source file written by a two-group fixed source simulation"""
    directory = tmp_path_factory.mktemp('mg_surface_source')
    plane = openmc.XPlane(0.0)
    model, _ = make_model('2g-thermal', directory, plane=plane)
    model.settings.run_mode = 'fixed source'
    model.settings.source = box_source([0.1, 1.0e6], [0.5, 0.5])
    model.settings.surf_source_write = {
        'surface_ids': [plane.id], 'max_particles': 100_000}
    run(model, cwd=directory)
    return directory / 'surface_source.h5'


def copy_source(src, dst, energy_mode=None, E=None, time=None):
    """Copy a source file, optionally changing its marker or its sites"""
    shutil.copy(src, dst)
    with h5py.File(dst, 'r+') as fh:
        if energy_mode is not None:
            fh.attrs['energy_mode'] = np.bytes_(energy_mode)
        if E is not None or time is not None:
            bank = fh['source_bank'][()]
            if E is not None:
                bank['E'] = E
            if time is not None:
                bank['time'] = time
            fh['source_bank'][...] = bank


def read_bank(path):
    with h5py.File(path, 'r') as fh:
        return fh['source_bank'][()]


def site_keys(sites):
    """Position, direction and weight of each site, as tuples"""
    if sites.dtype['r'].names:
        r = np.column_stack([sites['r'][c] for c in 'xyz'])
        u = np.column_stack([sites['u'][c] for c in 'xyz'])
    else:
        r, u = sites['r'], sites['u']
    return [tuple(x) + tuple(y) + (w,) for x, y, w in
            zip(r.tolist(), u.tolist(), sites['wgt'].tolist())]


def file_energies(path):
    """Map the position, direction and weight of each file site to its E"""
    bank = read_bank(path)
    return dict(zip(site_keys(bank), bank['E'].tolist()))


def write_ev_file(path):
    """Write 500 sites at 0.1 eV and 500 at 1 MeV, uniform in the cube"""
    rng = np.random.default_rng(1)
    positions = rng.uniform(-5.0, 5.0, (1000, 3))
    particles = [
        openmc.SourceParticle(r=tuple(r), u=(0.0, 0.0, 1.0),
                              E=0.1 if i < 500 else 1.0e6)
        for i, r in enumerate(positions)
    ]
    openmc.write_source_file(particles, path)


def count_bounds_warnings(text):
    return len(re.findall(r'WARNING: Energy bounds', ' '.join(text.split())))


# ---------------------------------------------------------------------------
# Tests that run OpenMC in a subprocess


@pytest.mark.parametrize('name', LIBRARIES)
def test_source_points_marked(name, mg_files):
    directory = mg_files(name)
    for filename in ('source.3.h5', 'source.h5'):
        with h5py.File(directory / filename, 'r') as fh:
            assert set(fh.attrs) == {'filetype', 'version', 'energy_mode'}
            assert fh.attrs['filetype'] == b'source'
            assert fh.attrs['energy_mode'] == b'multi-group'
            bank = fh['source_bank'][()]
        assert bank.size == 1000
        assert np.all(bank['E'] == birth_group(name))
    assert np.array_equal(read_bank(directory / 'source.3.h5'),
                          read_bank(directory / 'source.h5'))


@pytest.mark.parametrize('name', LIBRARIES)
def test_fixed_source_from_file(name, mg_files, run_in_tmpdir):
    shutil.copy(mg_files(name) / 'source.3.h5', 'source.3.h5')
    model, _ = make_model(name)
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.FileSource('source.3.h5')
    flux = run_flux(model, name)

    # The flux lies only in the energy bin of the stored group (bins are in
    # increasing energy, groups in decreasing energy)
    stored_bin = flux.size - 1 - birth_group(name)
    assert flux[stored_bin] > 0.0
    assert np.all(np.delete(flux, stored_bin) == 0.0)


def test_eigenvalue_from_file(mg_files, run_in_tmpdir):
    # The file holds sites in the low group, the only fissile one
    shutil.copy(mg_files('2g-thermal') / 'source.3.h5', 'source.3.h5')
    model, _ = make_model('2g-thermal')
    model.settings.run_mode = 'eigenvalue'
    model.settings.inactive = 1
    model.settings.source = openmc.FileSource('source.3.h5')
    flux = run_flux(model, '2g-thermal')
    assert flux[0] > 0.0
    assert flux[1] == 0.0


@pytest.mark.parametrize('energy_mode', [None, 'continuous-energy'])
def test_ev_file_transport(energy_mode, run_in_tmpdir):
    write_ev_file('ev.h5')
    with h5py.File('ev.h5', 'r') as fh:
        assert set(fh.attrs) == {'filetype', 'version'}
    copy_source('ev.h5', 'source.h5', energy_mode=energy_mode)
    model, _ = make_model('2g-thermal')
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.FileSource('source.h5')
    # Sites at 0.1 eV start in the low group and sites at 1 MeV in the high
    # group
    assert np.all(run_flux(model, '2g-thermal') > 0.0)


def test_unknown_energy_mode(mg_files, run_in_tmpdir):
    copy_source(mg_files('2g-thermal') / 'source.3.h5', 'source.h5',
                energy_mode='foo')
    model, _ = make_model('2g-thermal')
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.FileSource('source.h5')
    with pytest.raises(RuntimeError, match="'foo'"):
        run(model)


@pytest.mark.parametrize('E, pattern', [
    (2.0, r'group index 2(\.0)?\b'),
    (0.5, r'group index 0\.5\b'),
    (-1.0, r'group index -1(\.0)?\b'),
    (math.nan, r'group index nan\b'),
], ids=['G', 'fraction', 'negative', 'nan'])
def test_invalid_group_index(E, pattern, mg_files, run_in_tmpdir):
    copy_source(mg_files('2g-thermal') / 'source.3.h5', 'source.h5', E=E)
    model, _ = make_model('2g-thermal')
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.FileSource('source.h5')
    with pytest.raises(RuntimeError) as excinfo:
        run(model)
    message = str(excinfo.value)
    assert re.search(pattern, message, re.IGNORECASE)
    assert re.search(r'\b2 energy groups', message)


@pytest.mark.parametrize('filename', ['source.3.h5', 'statepoint.3.h5'])
def test_continuous_energy_refuses(filename, mg_files, run_in_tmpdir):
    statepoint = filename.startswith('statepoint')
    shutil.copy(mg_files('2g-thermal', statepoint) / filename, filename)

    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    mat.set_density('g/cm3', 0.1)
    cube = openmc.model.RectangularParallelepiped(
        -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, boundary_type='vacuum')
    model = openmc.Model()
    model.materials = openmc.Materials([mat])
    model.geometry = openmc.Geometry([openmc.Cell(fill=mat, region=-cube)])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 100
    model.settings.batches = 2
    model.settings.source = openmc.FileSource(filename)
    with pytest.raises(RuntimeError,
                       match='multigroup energy group indices'):
        run(model)


def test_continuous_energy_source_points_unmarked(run_in_tmpdir):
    fuel = openmc.Material()
    fuel.add_nuclide('U235', 1.0)
    fuel.add_nuclide('H1', 100.0)
    fuel.set_density('g/cm3', 1.0)
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    model = openmc.Model()
    model.materials = openmc.Materials([fuel])
    model.geometry = openmc.Geometry([openmc.Cell(fill=fuel, region=-sphere)])
    model.settings.particles = 500
    model.settings.batches = 3
    model.settings.inactive = 1
    model.settings.sourcepoint = {
        'batches': [3], 'separate': True, 'overwrite': True}
    run(model)
    for filename in ('source.3.h5', 'source.h5'):
        with h5py.File(filename, 'r') as fh:
            assert set(fh.attrs) == {'filetype', 'version'}


def test_surface_source_file_unmarked(surface_source_file):
    with h5py.File(surface_source_file, 'r') as fh:
        assert set(fh.attrs) == {'filetype', 'version'}
        E = fh['source_bank']['E']
    # Surface source sites hold the average energy of their group in [eV]
    assert set(np.unique(E).tolist()) == {LOW_AVERAGE, HIGH_AVERAGE}


def test_surface_source_transport(surface_source_file, run_in_tmpdir):
    shutil.copy(surface_source_file, 'surface_source.h5')
    model, _ = make_model('2g-thermal', plane=openmc.XPlane(0.0))
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.FileSource('surface_source.h5')
    assert np.all(run_flux(model, '2g-thermal') > 0.0)


def test_source_points_with_mpi(run_in_tmpdir, capsys):
    # With --mpi, the writer and the readers run on several processes
    Path('writer').mkdir()
    write_source_points('writer', '2g-thermal', mpi=True)
    with h5py.File('writer/source.3.h5', 'r') as fh:
        assert fh.attrs['energy_mode'] == b'multi-group'
        bank = fh['source_bank'][()]
    assert bank.size == 1000
    assert np.all(bank['E'] == 1.0)

    for directory in ('reader', 'wide'):
        Path(directory).mkdir()
        shutil.copy('writer/source.3.h5', directory)

    model, _ = make_model('2g-thermal', 'reader')
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.FileSource('source.3.h5')
    flux = run_flux(model, '2g-thermal', cwd='reader', mpi=True)
    assert flux[0] > 0.0
    assert flux[1] == 0.0

    # Energy bounds are not applied to group indices, which one warning says
    model, _ = make_model('2g-thermal', 'wide')
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.FileSource(
        'source.3.h5', constraints={'energy_bounds': WIDE_BOUNDS})
    capsys.readouterr()
    run(model, cwd='wide', mpi=True, output=True)
    output = capsys.readouterr().out
    assert count_bounds_warnings(output) == 1
    if config['mpi']:
        assert re.search(rf'MPI Processes \| {config["mpi_np"]}\b', output)


# ---------------------------------------------------------------------------
# Tests that sample the source through openmc.lib


def sample(model, n_samples=1000):
    model.export_to_model_xml()
    with openmc.lib.run_in_memory(output=False):
        return openmc.lib.sample_external_source(
            n_samples, prn_seed=1, as_array=True)


@pytest.mark.parametrize('bounds', [None, WIDE_BOUNDS],
                         ids=['default-bounds', 'wide-bounds'])
@pytest.mark.parametrize('name', LIBRARIES)
def test_sample_file_source(name, bounds, mg_files, run_in_tmpdir):
    shutil.copy(mg_files(name) / 'source.3.h5', 'source.3.h5')
    model, _ = make_model(name)
    constraints = {'energy_bounds': bounds} if bounds else None
    model.settings.source = openmc.FileSource(
        'source.3.h5', constraints=constraints)
    sites = sample(model)

    # Every sampled site is a file site and keeps its stored group
    stored = file_energies('source.3.h5')
    assert all(key in stored for key in site_keys(sites))
    assert np.all(sites['E'] == birth_group(name))


@pytest.mark.parametrize('bounds, n_warnings', [
    (None, 0),
    ((0.0, sys.float_info.max), 0),
    (WIDE_BOUNDS, 1),
    ((0.0, 1.0e7), 1),
    ((-1.0, sys.float_info.max), 1),
], ids=['none', 'default', 'wide', 'upper', 'lower'])
def test_energy_bounds_warning(bounds, n_warnings, mg_files, run_in_tmpdir,
                               capfd):
    shutil.copy(mg_files('2g-thermal') / 'source.3.h5', 'source.3.h5')
    model, _ = make_model('2g-thermal')
    constraints = {'energy_bounds': bounds} if bounds else None
    model.settings.source = openmc.FileSource(
        'source.3.h5', constraints=constraints)
    model.export_to_model_xml()
    capfd.readouterr()
    with openmc.lib.run_in_memory(output=False):
        pass
    assert count_bounds_warnings(capfd.readouterr().err) == n_warnings


def test_domain_constraint(mg_files, run_in_tmpdir):
    shutil.copy(mg_files('2g-thermal') / 'source.3.h5', 'source.3.h5')
    model, (inner, _) = make_model('2g-thermal')
    model.settings.source_rejection_fraction = 0.001
    model.settings.source = openmc.FileSource(
        'source.3.h5', constraints={'domains': [inner]})
    sites = sample(model)
    assert np.all(np.abs(sites['r']) < 2.0)
    assert np.all(sites['E'] == 1.0)


def test_time_constraint(mg_files, run_in_tmpdir):
    # Time bounds still apply to sites that hold group indices
    time = np.where(np.arange(1000) % 2 == 0, 1.0, 3.0)
    copy_source(mg_files('2g-thermal') / 'source.3.h5', 'source.h5',
                time=time)
    model, _ = make_model('2g-thermal')
    model.settings.source = openmc.FileSource(
        'source.h5', constraints={'time_bounds': [0.0, 2.0]})
    sites = sample(model)
    assert np.all(sites['time'] == 1.0)
    assert np.all(sites['E'] == 1.0)


def test_statepoint_source(mg_files, run_in_tmpdir):
    path = mg_files('2g-thermal', statepoint=True) / 'statepoint.3.h5'
    shutil.copy(path, 'statepoint.3.h5')
    with h5py.File('statepoint.3.h5', 'r') as fh:
        assert fh['energy_mode'][()] == b'multi-group'
        assert np.all(fh['source_bank']['E'] == 1.0)
    model, _ = make_model('2g-thermal')
    model.settings.source = openmc.FileSource('statepoint.3.h5')
    sites = sample(model)
    assert np.all(sites['E'] == 1.0)


@pytest.mark.parametrize('energy_mode', [None, 'continuous-energy'])
def test_ev_file_sampling(energy_mode, run_in_tmpdir):
    write_ev_file('ev.h5')
    copy_source('ev.h5', 'source.h5', energy_mode=energy_mode)
    model, _ = make_model('2g-thermal')
    model.settings.source = openmc.FileSource('source.h5')
    sites = sample(model)

    # The energies in eV are converted to groups: 0.1 eV is in group 1
    stored = file_energies('source.h5')
    expected = [1.0 if stored[key] == 0.1 else 0.0 for key in site_keys(sites)]
    assert np.array_equal(sites['E'], expected)
    assert 0 < np.count_nonzero(sites['E']) < sites.size


def test_ev_file_energy_bounds(run_in_tmpdir, capfd):
    # Energy bounds still apply to files that hold energies, with no warning
    write_ev_file('source.h5')
    model, _ = make_model('2g-thermal')
    model.settings.source = openmc.FileSource(
        'source.h5', constraints={'energy_bounds': [1.0, 3.0e7]})
    capfd.readouterr()
    sites = sample(model)
    assert np.all(sites['E'] == 0.0)
    assert count_bounds_warnings(capfd.readouterr().err) == 0


def test_surface_source_sampling(surface_source_file, run_in_tmpdir):
    shutil.copy(surface_source_file, 'surface_source.h5')
    model, _ = make_model('2g-thermal', plane=openmc.XPlane(0.0))
    model.settings.source = openmc.FileSource('surface_source.h5')
    sites = sample(model)

    # The group averages in eV map back to their groups
    stored = file_energies('surface_source.h5')
    expected = [1.0 if stored[key] == LOW_AVERAGE else 0.0
                for key in site_keys(sites)]
    assert np.array_equal(sites['E'], expected)
    assert 0 < np.count_nonzero(sites['E']) < sites.size


def test_mixed_sources(mg_files, run_in_tmpdir):
    # A file of group indices does not change how another source is read
    shutil.copy(mg_files('2g-thermal') / 'source.3.h5', 'source.3.h5')
    model, _ = make_model('2g-thermal')
    point = (4.5, 4.5, 4.5)
    model.settings.source = [
        openmc.FileSource('source.3.h5'),
        openmc.IndependentSource(
            space=openmc.stats.Point(point),
            energy=openmc.stats.Discrete([1.0e6], [1.0])),
    ]
    sites = sample(model)

    from_point = np.all(sites['r'] == point, axis=1)
    assert 0 < np.count_nonzero(from_point) < sites.size
    assert np.all(sites['E'][from_point] == 0.0)
    stored = file_energies('source.3.h5')
    file_sites = sites[~from_point]
    assert all(key in stored for key in site_keys(file_sites))
    assert np.all(file_sites['E'] == 1.0)
