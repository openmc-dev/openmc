from collections.abc import Mapping
import os

import numpy as np
import pytest
import openmc
import openmc.exceptions as exc
import openmc.lib

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

    # Add an energy function tally
    energyfunc_tally = openmc.Tally()
    energyfunc_filter = openmc.EnergyFunctionFilter(
        [0.0, 20e6], [0.0, 20e6])
    energyfunc_tally.scores = ['fission']
    energyfunc_tally.filters = [energyfunc_filter]
    pincell.tallies.append(energyfunc_tally)

    # Write XML files in tmpdir
    with cdtemp():
        pincell.export_to_xml()
        yield


@pytest.fixture(scope='module')
def lib_init(pincell_model, mpi_intracomm):
    openmc.lib.init(intracomm=mpi_intracomm)
    yield
    openmc.lib.finalize()


@pytest.fixture(scope='module')
def lib_simulation_init(lib_init):
    openmc.lib.simulation_init()
    yield


@pytest.fixture(scope='module')
def lib_run(lib_simulation_init):
    openmc.lib.run()


def test_cell_mapping(lib_init):
    cells = openmc.lib.cells
    assert isinstance(cells, Mapping)
    assert len(cells) == 3
    for cell_id, cell in cells.items():
        assert isinstance(cell, openmc.lib.Cell)
        assert cell_id == cell.id


def test_cell(lib_init):
    cell = openmc.lib.cells[1]
    assert isinstance(cell.fill, openmc.lib.Material)
    cell.fill = openmc.lib.materials[1]
    assert str(cell) == 'Cell[0]'
    assert cell.name == "Fuel"
    cell.name = "Not fuel"
    assert cell.name == "Not fuel"

def test_cell_temperature(lib_init):
    cell = openmc.lib.cells[1]
    cell.set_temperature(100.0, 0)
    assert cell.get_temperature(0) == 100.0
    cell.set_temperature(200)
    assert cell.get_temperature() == 200.0


def test_new_cell(lib_init):
    with pytest.raises(exc.AllocationError):
        openmc.lib.Cell(1)
    new_cell = openmc.lib.Cell()
    new_cell_with_id = openmc.lib.Cell(10)
    assert len(openmc.lib.cells) == 5


def test_material_mapping(lib_init):
    mats = openmc.lib.materials
    assert isinstance(mats, Mapping)
    assert len(mats) == 3
    for mat_id, mat in mats.items():
        assert isinstance(mat, openmc.lib.Material)
        assert mat_id == mat.id


def test_material(lib_init):
    m = openmc.lib.materials[3]
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
    assert m.name == "Hot borated water"
    m.name = "Not hot borated water"
    assert m.name == "Not hot borated water"

def test_material_add_nuclide(lib_init):
    m = openmc.lib.materials[3]
    m.add_nuclide('Xe135', 1e-12)
    assert m.nuclides[-1] == 'Xe135'
    assert m.densities[-1] == 1e-12


def test_new_material(lib_init):
    with pytest.raises(exc.AllocationError):
        openmc.lib.Material(1)
    new_mat = openmc.lib.Material()
    new_mat_with_id = openmc.lib.Material(10)
    assert len(openmc.lib.materials) == 5


def test_nuclide_mapping(lib_init):
    nucs = openmc.lib.nuclides
    assert isinstance(nucs, Mapping)
    assert len(nucs) == 13
    for name, nuc in nucs.items():
        assert isinstance(nuc, openmc.lib.Nuclide)
        assert name == nuc.name


def test_settings(lib_init):
    settings = openmc.lib.settings
    assert settings.batches == 10
    settings.batches = 10
    assert settings.inactive == 5
    assert settings.generations_per_batch == 1
    assert settings.particles == 100
    assert settings.seed == 1
    settings.seed = 11


def test_tally_mapping(lib_init):
    tallies = openmc.lib.tallies
    assert isinstance(tallies, Mapping)
    assert len(tallies) == 3
    for tally_id, tally in tallies.items():
        assert isinstance(tally, openmc.lib.Tally)
        assert tally_id == tally.id


def test_energy_function_filter(lib_init):
    """Test special __new__ and __init__ for EnergyFunctionFilter"""
    efunc = openmc.lib.EnergyFunctionFilter([0.0, 1.0], [0.0, 2.0])
    assert len(efunc.energy) == 2
    assert (efunc.energy == [0.0, 1.0]).all()
    assert len(efunc.y) == 2
    assert (efunc.y == [0.0, 2.0]).all()


def test_tally(lib_init):
    t = openmc.lib.tallies[1]
    assert t.type == 'volume'
    assert len(t.filters) == 2
    assert isinstance(t.filters[0], openmc.lib.MaterialFilter)
    assert isinstance(t.filters[1], openmc.lib.EnergyFilter)

    # Create new filter and replace existing
    with pytest.raises(exc.AllocationError):
        openmc.lib.MaterialFilter(uid=1)
    mats = openmc.lib.materials
    f = openmc.lib.MaterialFilter([mats[2], mats[1]])
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

    t2 = openmc.lib.tallies[2]
    assert len(t2.filters) == 2
    assert isinstance(t2.filters[0], openmc.lib.ZernikeFilter)
    assert isinstance(t2.filters[1], openmc.lib.CellFilter)
    assert len(t2.filters[1].bins) == 3
    assert t2.filters[0].order == 5

    t3 = openmc.lib.tallies[3]
    assert len(t3.filters) == 1
    t3_f = t3.filters[0]
    assert isinstance(t3_f, openmc.lib.EnergyFunctionFilter)
    assert len(t3_f.energy) == 2
    assert len(t3_f.y) == 2
    t3_f.set_data([0.0, 1.0, 2.0], [0.0, 1.0, 4.0])
    assert len(t3_f.energy) == 3
    assert len(t3_f.y) == 3


def test_new_tally(lib_init):
    with pytest.raises(exc.AllocationError):
        openmc.lib.Material(1)
    new_tally = openmc.lib.Tally()
    new_tally.scores = ['flux']
    new_tally_with_id = openmc.lib.Tally(10)
    new_tally_with_id.scores = ['flux']
    assert len(openmc.lib.tallies) == 5


def test_tally_activate(lib_simulation_init):
    t = openmc.lib.tallies[1]
    assert not t.active
    t.active = True
    assert t.active


def test_tally_writable(lib_simulation_init):
    t = openmc.lib.tallies[1]
    assert t.writable
    t.writable = False
    assert not t.writable
    # Revert tally to writable state for lib_run fixtures
    t.writable = True


def test_tally_results(lib_run):
    t = openmc.lib.tallies[1]
    assert t.num_realizations == 10  # t was made active in test_tally_active
    assert np.all(t.mean >= 0)
    nonzero = (t.mean > 0.0)
    assert np.all(t.std_dev[nonzero] >= 0)
    assert np.all(t.ci_width()[nonzero] >= 1.95*t.std_dev[nonzero])

    t2 = openmc.lib.tallies[2]
    n = 5
    assert t2.mean.size == (n + 1) * (n + 2) // 2 * 3 # Number of Zernike coeffs * 3 cells


def test_global_tallies(lib_run):
    assert openmc.lib.num_realizations() == 5
    gt = openmc.lib.global_tallies()
    for mean, std_dev in gt:
        assert mean >= 0


def test_statepoint(lib_run):
    openmc.lib.statepoint_write('test_sp.h5')
    assert os.path.exists('test_sp.h5')


def test_source_bank(lib_run):
    source = openmc.lib.source_bank()
    assert np.all(source['E'] > 0.0)
    assert np.all(source['wgt'] == 1.0)
    assert np.allclose(np.linalg.norm(source['u'], axis=1), 1.0)


def test_by_batch(lib_run):
    openmc.lib.hard_reset()

    # Running next batch before simulation is initialized should raise an
    # exception
    with pytest.raises(exc.AllocationError):
        openmc.lib.next_batch()

    openmc.lib.simulation_init()
    try:
        for _ in openmc.lib.iter_batches():
            # Make sure we can get k-effective during inactive/active batches
            mean, std_dev = openmc.lib.keff()
            assert 0.0 < mean < 2.5
            assert std_dev > 0.0
        assert openmc.lib.num_realizations() == 5

        for i in range(3):
            openmc.lib.next_batch()
        assert openmc.lib.num_realizations() == 8

    finally:
        openmc.lib.simulation_finalize()


def test_reset(lib_run):
    # Init and run 10 batches.
    openmc.lib.hard_reset()
    openmc.lib.simulation_init()
    try:
        for i in range(10):
            openmc.lib.next_batch()

        # Make sure there are 5 realizations for the 5 active batches.
        assert openmc.lib.num_realizations() == 5
        assert openmc.lib.tallies[2].num_realizations == 5
        _, keff_sd1 = openmc.lib.keff()
        tally_sd1 = openmc.lib.tallies[2].std_dev[0]

        # Reset and run 3 more batches.  Check the number of realizations.
        openmc.lib.reset()
        for i in range(3):
            openmc.lib.next_batch()
        assert openmc.lib.num_realizations() == 3
        assert openmc.lib.tallies[2].num_realizations == 3

        # Check the tally std devs to make sure results were cleared.
        _, keff_sd2 = openmc.lib.keff()
        tally_sd2 = openmc.lib.tallies[2].std_dev[0]
        assert keff_sd2 > keff_sd1
        assert tally_sd2 > tally_sd1

    finally:
        openmc.lib.simulation_finalize()


def test_reproduce_keff(lib_init):
    # Get k-effective after run
    openmc.lib.hard_reset()
    openmc.lib.run()
    keff0 = openmc.lib.keff()

    # Reset, run again, and get k-effective again. they should match
    openmc.lib.hard_reset()
    openmc.lib.run()
    keff1 = openmc.lib.keff()
    assert keff0 == pytest.approx(keff1)


def test_find_cell(lib_init):
    cell, instance = openmc.lib.find_cell((0., 0., 0.))
    assert cell is openmc.lib.cells[1]
    cell, instance = openmc.lib.find_cell((0.4, 0., 0.))
    assert cell is openmc.lib.cells[2]
    with pytest.raises(exc.GeometryError):
        openmc.lib.find_cell((100., 100., 100.))


def test_find_material(lib_init):
    mat = openmc.lib.find_material((0., 0., 0.))
    assert mat is openmc.lib.materials[1]
    mat = openmc.lib.find_material((0.4, 0., 0.))
    assert mat is openmc.lib.materials[2]


def test_mesh(lib_init):
    mesh = openmc.lib.RegularMesh()
    mesh.dimension = (2, 3, 4)
    assert mesh.dimension == (2, 3, 4)
    with pytest.raises(exc.AllocationError):
        mesh2 = openmc.lib.RegularMesh(mesh.id)

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

    meshes = openmc.lib.meshes
    assert isinstance(meshes, Mapping)
    assert len(meshes) == 1
    for mesh_id, mesh in meshes.items():
        assert isinstance(mesh, openmc.lib.RegularMesh)
        assert mesh_id == mesh.id

    mf = openmc.lib.MeshFilter(mesh)
    assert mf.mesh == mesh

    msf = openmc.lib.MeshSurfaceFilter(mesh)
    assert msf.mesh == mesh


def test_restart(lib_init, mpi_intracomm):
    # Finalize and re-init to make internal state consistent with XML.
    openmc.lib.hard_reset()
    openmc.lib.finalize()
    openmc.lib.init(intracomm=mpi_intracomm)
    openmc.lib.simulation_init()

    # Run for 7 batches then write a statepoint.
    for i in range(7):
        openmc.lib.next_batch()
    openmc.lib.statepoint_write('restart_test.h5', True)

    # Run 3 more batches and copy the keff.
    for i in range(3):
        openmc.lib.next_batch()
    keff0 = openmc.lib.keff()

    # Restart the simulation from the statepoint and the 3 remaining active batches.
    openmc.lib.simulation_finalize()
    openmc.lib.hard_reset()
    openmc.lib.finalize()
    openmc.lib.init(args=('-r', 'restart_test.h5'))
    openmc.lib.simulation_init()
    for i in range(3):
        openmc.lib.next_batch()
    keff1 = openmc.lib.keff()
    openmc.lib.simulation_finalize()

    # Compare the keff values.
    assert keff0 == pytest.approx(keff1)


def test_load_nuclide(lib_init):
    # load multiple nuclides
    openmc.lib.load_nuclide('H3')
    assert 'H3' in openmc.lib.nuclides
    openmc.lib.load_nuclide('Pu239')
    assert 'Pu239' in openmc.lib.nuclides
    # load non-existent nuclide
    with pytest.raises(exc.DataError):
        openmc.lib.load_nuclide('Pu3')


def test_id_map(lib_init):
    expected_ids = np.array([[(3, 3), (2, 2), (3, 3)],
                             [(2, 2), (1, 1), (2, 2)],
                             [(3, 3), (2, 2), (3, 3)]], dtype='int32')

    # create a plot object
    s = openmc.lib.plot._PlotBase()
    s.width = 1.26
    s.height = 1.26
    s.v_res = 3
    s.h_res = 3
    s.origin = (0.0, 0.0, 0.0)
    s.basis = 'xy'
    s.level = -1

    ids = openmc.lib.plot.id_map(s)
    assert np.array_equal(expected_ids, ids)

def test_property_map(lib_init):
    expected_properties = np.array(
        [[(293.6, 0.740582), (293.6, 6.55), (293.6, 0.740582)],
         [ (293.6, 6.55), (293.6, 10.29769),  (293.6, 6.55)],
         [(293.6, 0.740582), (293.6, 6.55), (293.6, 0.740582)]], dtype='float')

    # create a plot object
    s = openmc.lib.plot._PlotBase()
    s.width = 1.26
    s.height = 1.26
    s.v_res = 3
    s.h_res = 3
    s.origin = (0.0, 0.0, 0.0)
    s.basis = 'xy'
    s.level = -1

    properties = openmc.lib.plot.property_map(s)
    assert np.allclose(expected_properties, properties, atol=1e-04)


def test_position(lib_init):

    pos = openmc.lib.plot._Position(1.0, 2.0, 3.0)

    assert tuple(pos) == (1.0, 2.0, 3.0)

    pos[0] = 1.3
    pos[1] = 2.3
    pos[2] = 3.3

    assert tuple(pos) == (1.3, 2.3, 3.3)


def test_global_bounding_box(lib_init):
    expected_llc = (-0.63, -0.63, -np.inf)
    expected_urc = (0.63, 0.63, np.inf)

    llc, urc = openmc.lib.global_bounding_box()

    assert tuple(llc) == expected_llc
    assert tuple(urc) == expected_urc
