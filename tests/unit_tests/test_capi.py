#!/usr/bin/env python

from collections import Mapping
import os

import numpy as np
import pytest
import openmc
import openmc.capi


@pytest.fixture(scope='module')
def pincell_model():
    """Set up a model to test with and delete files when done"""
    pincell = openmc.examples.pwr_pin_cell()

    # Add a tally
    filter1 = openmc.MaterialFilter(pincell.materials)
    filter2 = openmc.EnergyFilter([0.0, 1.0, 1.0e3, 20.0e6])
    mat_tally = openmc.Tally()
    mat_tally.filters = [filter1, filter2]
    mat_tally.nuclides = ['U235', 'U238']
    mat_tally.scores = ['total', 'elastic', '(n,gamma)']
    pincell.tallies.append(mat_tally)

    # Write XML files
    pincell.export_to_xml()

    yield

    # Delete generated files
    files = ['geometry.xml', 'materials.xml', 'settings.xml', 'tallies.xml',
             'statepoint.10.h5', 'summary.h5', 'test_sp.h5']
    for f in files:
        if os.path.exists(f):
            os.remove(f)


@pytest.fixture(scope='module')
def capi_init(pincell_model):
    openmc.capi.init()
    yield
    openmc.capi.finalize()


@pytest.fixture(scope='module')
def capi_run(capi_init):
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
    assert str(cell) == 'Cell[1]'


def test_new_cell(capi_init):
    with pytest.raises(openmc.capi.AllocationError):
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

    rho = 2.25e-2
    m.set_density(rho)
    assert sum(m.densities) == pytest.approx(rho)


def test_new_material(capi_init):
    with pytest.raises(openmc.capi.AllocationError):
        openmc.capi.Material(1)
    new_mat = openmc.capi.Material()
    new_mat_with_id = openmc.capi.Material(10)
    assert len(openmc.capi.materials) == 5


def test_nuclide_mapping(capi_init):
    nucs = openmc.capi.nuclides
    assert isinstance(nucs, Mapping)
    assert len(nucs) == 12
    for name, nuc in nucs.items():
        assert isinstance(nuc, openmc.capi.Nuclide)
        assert name == nuc.name


def test_load_nuclide(capi_init):
    openmc.capi.load_nuclide('Pu239')
    with pytest.raises(openmc.capi.DataError):
        openmc.capi.load_nuclide('Pu3')


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
    assert len(tallies) == 1
    for tally_id, tally in tallies.items():
        assert isinstance(tally, openmc.capi.Tally)
        assert tally_id == tally.id


def test_tally(capi_init):
    t = openmc.capi.tallies[1]
    t.id = 1
    assert len(t.filters) == 2
    assert isinstance(t.filters[0], openmc.capi.MaterialFilter)
    assert isinstance(t.filters[1], openmc.capi.EnergyFilter)

    # Create new filter and replace existing
    with pytest.raises(openmc.capi.AllocationError):
        openmc.capi.MaterialFilter(uid=1)
    mats = openmc.capi.materials
    f = openmc.capi.MaterialFilter([mats[2], mats[1]])
    t.filters = [f]
    assert t.filters == [f]

    assert t.nuclides == ['U235', 'U238']
    with pytest.raises(openmc.capi.DataError):
        t.nuclides = ['Zr2']
    t.nuclides = ['U234', 'Zr90']
    assert t.nuclides == ['U234', 'Zr90']

    assert t.scores == ['total', '(n,elastic)', '(n,gamma)']
    new_scores = ['scatter', 'fission', 'nu-fission', '(n,2n)']
    t.scores = new_scores
    assert t.scores == new_scores


def test_new_tally(capi_init):
    with pytest.raises(openmc.capi.AllocationError):
        openmc.capi.Material(1)
    new_tally = openmc.capi.Tally()
    new_tally.scores = ['flux']
    new_tally_with_id = openmc.capi.Tally(10)
    new_tally_with_id.scores = ['flux']
    assert len(openmc.capi.tallies) == 3


def test_tally_results(capi_run):
    t = openmc.capi.tallies[1]
    assert t.num_realizations == 5
    assert np.all(t.mean >= 0)
    nonzero = (t.mean > 0.0)
    assert np.all(t.std_dev[nonzero] >= 0)
    assert np.all(t.ci_width()[nonzero] >= 1.95*t.std_dev[nonzero])


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


def test_by_batch(capi_run):
    openmc.capi.hard_reset()

    # Running next batch before simulation is initialized should raise an
    # exception
    with pytest.raises(openmc.capi.AllocationError):
        openmc.capi.next_batch()

    openmc.capi.simulation_init()
    for _ in openmc.capi.iter_batches():
        # Make sure we can get k-effective during inactive/active batches
        mean, std_dev = openmc.capi.keff()
        assert 0.0 < mean < 2.5
        assert std_dev > 0.0
    assert openmc.capi.num_realizations() == 5

    for i in range(3):
        openmc.capi.next_batch()
    assert openmc.capi.num_realizations() == 8
    openmc.capi.simulation_finalize()


def test_find_cell(capi_init):
    cell, instance = openmc.capi.find_cell((0., 0., 0.))
    assert cell is openmc.capi.cells[1]
    cell, instance = openmc.capi.find_cell((0.4, 0., 0.))
    assert cell is openmc.capi.cells[2]
    with pytest.raises(openmc.capi.GeometryError):
        openmc.capi.find_cell((100., 100., 100.))


def test_find_material(capi_init):
    mat = openmc.capi.find_material((0., 0., 0.))
    assert mat is openmc.capi.materials[1]
    mat = openmc.capi.find_material((0.4, 0., 0.))
    assert mat is openmc.capi.materials[2]
