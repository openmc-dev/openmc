from collections.abc import Mapping
import os

import numpy as np
import pytest
import openmc
import openmc.exceptions as exc
import openmc.capi

from tests import cdtemp


@pytest.fixture(scope='module')
def pincell_model():
    """Set up a model to test with and delete files when done"""
    openmc.reset_auto_ids()
    pincell = openmc.examples.pwr_pin_cell()
    pincell.settings.verbosity = 1

    # Add a tally
    filter1 = openmc.MaterialFilter(pincell.materials)
    filter2 = openmc.EnergyFilter([0.0, 1.0, 1.0e3, 20.0e6])
    mat_tally = openmc.Tally()
    mat_tally.filters = [filter1, filter2]
    mat_tally.nuclides = ['U235', 'U238']
    mat_tally.scores = ['total', 'elastic', '(n,gamma)']
    pincell.tallies.append(mat_tally)

    # Add an expansion tally
    zernike_tally = openmc.Tally()
    filter3 = openmc.ZernikeFilter(5, r=.63)
    cells = pincell.geometry.root_universe.cells
    filter4 = openmc.CellFilter(list(cells.values()))
    zernike_tally.filters = [filter3, filter4]
    zernike_tally.scores = ['fission']
    pincell.tallies.append(zernike_tally)

    # Write XML files in tmpdir
    with cdtemp():
        pincell.export_to_xml()
        yield


@pytest.fixture(scope='module')
def capi_init(pincell_model):
    openmc.capi.init()
    yield
    openmc.capi.finalize()


@pytest.fixture(scope='module')
def capi_simulation_init(capi_init):
    openmc.capi.simulation_init()
    yield


@pytest.fixture(scope='module')
def capi_run(capi_simulation_init):
    openmc.capi.run()


def test_cell_mapping(capi_init):
    cells = openmc.capi.cells
    assert isinstance(cells, Mapping)
    assert len(cells) == 3
    for cell_id, cell in cells.items():
        assert isinstance(cell, openmc.capi.Cell)
        assert cell_id == cell.id


def test_cell(capi_init):
    cell = openmc.capi.cells[1]
    assert isinstance(cell.fill, openmc.capi.Material)
    cell.fill = openmc.capi.materials[1]
    assert str(cell) == 'Cell[0]'


def test_cell_temperature(capi_init):
    cell = openmc.capi.cells[1]
    cell.set_temperature(100.0, 0)
    assert cell.get_temperature(0) == 100.0
    cell.set_temperature(200)
    assert cell.get_temperature() == 200.0


def test_new_cell(capi_init):
    with pytest.raises(exc.AllocationError):
        openmc.capi.Cell(1)
    new_cell = openmc.capi.Cell()
    new_cell_with_id = openmc.capi.Cell(10)
    assert len(openmc.capi.cells) == 5


def test_material_mapping(capi_init):
    mats = openmc.capi.materials
    assert isinstance(mats, Mapping)
    assert len(mats) == 3
    for mat_id, mat in mats.items():
        assert isinstance(mat, openmc.capi.Material)
        assert mat_id == mat.id


def test_material(capi_init):
    m = openmc.capi.materials[3]
    assert m.nuclides == ['H1', 'O16', 'B10', 'B11']

    old_dens = m.densities
    test_dens = [1.0e-1, 2.0e-1, 2.5e-1, 1.0e-3]
    m.set_densities(m.nuclides, test_dens)
    assert m.densities == pytest.approx(test_dens)

    assert m.volume is None
    m.volume = 10.0
    assert m.volume == 10.0

    with pytest.raises(exc.OpenMCError):
        m.set_density(1.0, 'goblins')

    rho = 2.25e-2
    m.set_density(rho)
    assert sum(m.densities) == pytest.approx(rho)

    m.set_density(0.1, 'g/cm3')
    assert m.density == pytest.approx(0.1)


def test_material_add_nuclide(capi_init):
    m = openmc.capi.materials[3]
    m.add_nuclide('Xe135', 1e-12)
    assert m.nuclides[-1] == 'Xe135'
    assert m.densities[-1] == 1e-12


def test_new_material(capi_init):
    with pytest.raises(exc.AllocationError):
        openmc.capi.Material(1)
    new_mat = openmc.capi.Material()
    new_mat_with_id = openmc.capi.Material(10)
    assert len(openmc.capi.materials) == 5


def test_nuclide_mapping(capi_init):
    nucs = openmc.capi.nuclides
    assert isinstance(nucs, Mapping)
    assert len(nucs) == 13
    for name, nuc in nucs.items():
        assert isinstance(nuc, openmc.capi.Nuclide)
        assert name == nuc.name


def test_settings(capi_init):
    settings = openmc.capi.settings
    assert settings.batches == 10
    settings.batches = 10
    assert settings.inactive == 5
    assert settings.generations_per_batch == 1
    assert settings.particles == 100
    assert settings.seed == 1
    settings.seed = 11

    assert settings.run_mode == 'eigenvalue'
    settings.run_mode = 'volume'
    settings.run_mode = 'eigenvalue'


def test_tally_mapping(capi_init):
    tallies = openmc.capi.tallies
    assert isinstance(tallies, Mapping)
    assert len(tallies) == 2
    for tally_id, tally in tallies.items():
        assert isinstance(tally, openmc.capi.Tally)
        assert tally_id == tally.id


def test_tally(capi_init):
    t = openmc.capi.tallies[1]
    assert t.type == 'volume'
    assert len(t.filters) == 2
    assert isinstance(t.filters[0], openmc.capi.MaterialFilter)
    assert isinstance(t.filters[1], openmc.capi.EnergyFilter)

    # Create new filter and replace existing
    with pytest.raises(exc.AllocationError):
        openmc.capi.MaterialFilter(uid=1)
    mats = openmc.capi.materials
    f = openmc.capi.MaterialFilter([mats[2], mats[1]])
    assert f.bins[0] == mats[2]
    assert f.bins[1] == mats[1]
    t.filters = [f]
    assert t.filters == [f]

    assert t.nuclides == ['U235', 'U238']
    with pytest.raises(exc.DataError):
        t.nuclides = ['Zr2']
    t.nuclides = ['U234', 'Zr90']
    assert t.nuclides == ['U234', 'Zr90']

    assert t.scores == ['total', '(n,elastic)', '(n,gamma)']
    new_scores = ['scatter', 'fission', 'nu-fission', '(n,2n)']
    t.scores = new_scores
    assert t.scores == new_scores

    t2 = openmc.capi.tallies[2]
    assert len(t2.filters) == 2
    assert isinstance(t2.filters[0], openmc.capi.ZernikeFilter)
    assert isinstance(t2.filters[1], openmc.capi.CellFilter)
    assert len(t2.filters[1].bins) == 3
    assert t2.filters[0].order == 5


def test_new_tally(capi_init):
    with pytest.raises(exc.AllocationError):
        openmc.capi.Material(1)
    new_tally = openmc.capi.Tally()
    new_tally.scores = ['flux']
    new_tally_with_id = openmc.capi.Tally(10)
    new_tally_with_id.scores = ['flux']
    assert len(openmc.capi.tallies) == 4


def test_tally_activate(capi_simulation_init):
    t = openmc.capi.tallies[1]
    assert not t.active
    t.active = True
    assert t.active


def test_tally_results(capi_run):
    t = openmc.capi.tallies[1]
    assert t.num_realizations == 10  # t was made active in test_tally
    assert np.all(t.mean >= 0)
    nonzero = (t.mean > 0.0)
    assert np.all(t.std_dev[nonzero] >= 0)
    assert np.all(t.ci_width()[nonzero] >= 1.95*t.std_dev[nonzero])

    t2 = openmc.capi.tallies[2]
    n = 5
    assert t2.mean.size == (n + 1) * (n + 2) // 2 * 3 # Number of Zernike coeffs * 3 cells


def test_global_tallies(capi_run):
    assert openmc.capi.num_realizations() == 5
    gt = openmc.capi.global_tallies()
    for mean, std_dev in gt:
        assert mean >= 0


def test_statepoint(capi_run):
    openmc.capi.statepoint_write('test_sp.h5')
    assert os.path.exists('test_sp.h5')


def test_source_bank(capi_run):
    source = openmc.capi.source_bank()
    assert np.all(source['E'] > 0.0)
    assert np.all(source['wgt'] == 1.0)
    assert np.allclose(np.linalg.norm(source['u'], axis=1), 1.0)


def test_by_batch(capi_run):
    openmc.capi.hard_reset()

    # Running next batch before simulation is initialized should raise an
    # exception
    with pytest.raises(exc.AllocationError):
        openmc.capi.next_batch()

    openmc.capi.simulation_init()
    try:
        for _ in openmc.capi.iter_batches():
            # Make sure we can get k-effective during inactive/active batches
            mean, std_dev = openmc.capi.keff()
            assert 0.0 < mean < 2.5
            assert std_dev > 0.0
        assert openmc.capi.num_realizations() == 5

        for i in range(3):
            openmc.capi.next_batch()
        assert openmc.capi.num_realizations() == 8

    finally:
        openmc.capi.simulation_finalize()


def test_reset(capi_run):
    # Init and run 10 batches.
    openmc.capi.hard_reset()
    openmc.capi.simulation_init()
    try:
        for i in range(10):
            openmc.capi.next_batch()

        # Make sure there are 5 realizations for the 5 active batches.
        assert openmc.capi.num_realizations() == 5
        assert openmc.capi.tallies[2].num_realizations == 5
        _, keff_sd1 = openmc.capi.keff()
        tally_sd1 = openmc.capi.tallies[2].std_dev[0]

        # Reset and run 3 more batches.  Check the number of realizations.
        openmc.capi.reset()
        for i in range(3):
            openmc.capi.next_batch()
        assert openmc.capi.num_realizations() == 3
        assert openmc.capi.tallies[2].num_realizations == 3

        # Check the tally std devs to make sure results were cleared.
        _, keff_sd2 = openmc.capi.keff()
        tally_sd2 = openmc.capi.tallies[2].std_dev[0]
        assert keff_sd2 > keff_sd1
        assert tally_sd2 > tally_sd1

    finally:
        openmc.capi.simulation_finalize()


def test_reproduce_keff(capi_init):
    # Get k-effective after run
    openmc.capi.hard_reset()
    openmc.capi.run()
    keff0 = openmc.capi.keff()

    # Reset, run again, and get k-effective again. they should match
    openmc.capi.hard_reset()
    openmc.capi.run()
    keff1 = openmc.capi.keff()
    assert keff0 == pytest.approx(keff1)


def test_find_cell(capi_init):
    cell, instance = openmc.capi.find_cell((0., 0., 0.))
    assert cell is openmc.capi.cells[1]
    cell, instance = openmc.capi.find_cell((0.4, 0., 0.))
    assert cell is openmc.capi.cells[2]
    with pytest.raises(exc.GeometryError):
        openmc.capi.find_cell((100., 100., 100.))


def test_find_material(capi_init):
    mat = openmc.capi.find_material((0., 0., 0.))
    assert mat is openmc.capi.materials[1]
    mat = openmc.capi.find_material((0.4, 0., 0.))
    assert mat is openmc.capi.materials[2]


def test_mesh(capi_init):
    mesh = openmc.capi.RegularMesh()
    mesh.dimension = (2, 3, 4)
    assert mesh.dimension == (2, 3, 4)
    with pytest.raises(exc.AllocationError):
        mesh2 = openmc.capi.RegularMesh(mesh.id)

    # Make sure each combination of parameters works
    ll = (0., 0., 0.)
    ur = (10., 10., 10.)
    width = (1., 1., 1.)
    mesh.set_parameters(lower_left=ll, upper_right=ur)
    assert mesh.lower_left == pytest.approx(ll)
    assert mesh.upper_right == pytest.approx(ur)
    mesh.set_parameters(lower_left=ll, width=width)
    assert mesh.lower_left == pytest.approx(ll)
    assert mesh.width == pytest.approx(width)
    mesh.set_parameters(upper_right=ur, width=width)
    assert mesh.upper_right == pytest.approx(ur)
    assert mesh.width == pytest.approx(width)

    meshes = openmc.capi.meshes
    assert isinstance(meshes, Mapping)
    assert len(meshes) == 1
    for mesh_id, mesh in meshes.items():
        assert isinstance(mesh, openmc.capi.RegularMesh)
        assert mesh_id == mesh.id

    mf = openmc.capi.MeshFilter(mesh)
    assert mf.mesh == mesh

    msf = openmc.capi.MeshSurfaceFilter(mesh)
    assert msf.mesh == mesh


def test_restart(capi_init):
    # Finalize and re-init to make internal state consistent with XML.
    openmc.capi.hard_reset()
    openmc.capi.finalize()
    openmc.capi.init()
    openmc.capi.simulation_init()

    # Run for 7 batches then write a statepoint.
    for i in range(7):
        openmc.capi.next_batch()
    openmc.capi.statepoint_write('restart_test.h5', True)

    # Run 3 more batches and copy the keff.
    for i in range(3):
        openmc.capi.next_batch()
    keff0 = openmc.capi.keff()

    # Restart the simulation from the statepoint and the 3 remaining active batches.
    openmc.capi.simulation_finalize()
    openmc.capi.hard_reset()
    openmc.capi.finalize()
    openmc.capi.init(args=('-r', 'restart_test.h5'))
    openmc.capi.simulation_init()
    for i in range(3):
        openmc.capi.next_batch()
    keff1 = openmc.capi.keff()
    openmc.capi.simulation_finalize()

    # Compare the keff values.
    assert keff0 == pytest.approx(keff1)


def test_load_nuclide(capi_init):
    # load multiple nuclides
    openmc.capi.load_nuclide('H3')
    assert 'H3' in openmc.capi.nuclides
    openmc.capi.load_nuclide('Pu239')
    assert 'Pu239' in openmc.capi.nuclides
    # load non-existent nuclide
    with pytest.raises(exc.DataError):
        openmc.capi.load_nuclide('Pu3')


def test_id_map(capi_init):
    expected_ids = np.array([[(3, 3), (2, 2), (3, 3)],
                             [(2, 2), (1, 1), (2, 2)],
                             [(3, 3), (2, 2), (3, 3)]], dtype='int32')

    # create a plot object
    s = openmc.capi.plot._PlotBase()
    s.width = 1.26
    s.height = 1.26
    s.v_res = 3
    s.h_res = 3
    s.origin = (0.0, 0.0, 0.0)
    s.basis = 'xy'
    s.level = -1

    ids = openmc.capi.plot.id_map(s)
    assert np.array_equal(expected_ids, ids)

def test_property_map(capi_init):
    expected_properties = np.array(
        [[(293.6, 0.740582), (293.6, 6.55), (293.6, 0.740582)],
         [ (293.6, 6.55), (293.6, 10.29769),  (293.6, 6.55)],
         [(293.6, 0.740582), (293.6, 6.55), (293.6, 0.740582)]], dtype='float')

    # create a plot object
    s = openmc.capi.plot._PlotBase()
    s.width = 1.26
    s.height = 1.26
    s.v_res = 3
    s.h_res = 3
    s.origin = (0.0, 0.0, 0.0)
    s.basis = 'xy'
    s.level = -1

    properties = openmc.capi.plot.property_map(s)
    assert np.allclose(expected_properties, properties, atol=1e-04)


def test_position(capi_init):

    pos = openmc.capi.plot._Position(1.0, 2.0, 3.0)

    assert tuple(pos) == (1.0, 2.0, 3.0)

    pos[0] = 1.3
    pos[1] = 2.3
    pos[2] = 3.3

    assert tuple(pos) == (1.3, 2.3, 3.3)
