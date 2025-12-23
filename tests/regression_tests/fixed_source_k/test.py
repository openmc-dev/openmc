from glob import glob
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config
import os
import glob

def create_universe():
    # Define materials
    heu = openmc.Material(name='HEU')
    heu.add_element('U', 1, enrichment=93.0, enrichment_type='wo')
    heu.set_density('g/cm3', 19.1)

    dep_uranium = openmc.Material(name='Depleted Uranium')
    dep_uranium.add_nuclide('U238', 1.0)
    dep_uranium.set_density('g/cm3', 19.1)

    mats = openmc.Materials([heu, dep_uranium])
    mats.export_to_xml()

    # Geometry
    fuel_radius = 4.8
    dep_uranium_radius = 12.2

    fuel_sphere = openmc.Sphere(r=fuel_radius)
    dep_uranium_sphere = openmc.Sphere(r=dep_uranium_radius, boundary_type='vacuum')

    fuel_region = -fuel_sphere
    dep_uranium_region = +fuel_sphere & -dep_uranium_sphere

    fuel_cell = openmc.Cell(name='Fuel', region=fuel_region)
    fuel_cell.fill = heu

    dep_uranium_cell = openmc.Cell(name='Depleted Uranium', region=dep_uranium_region)
    dep_uranium_cell.fill = dep_uranium

    universe = openmc.Universe(cells=[fuel_cell, dep_uranium_cell])
    return universe

def test_source():
    source_space = openmc.stats.Point((0,0,0))
    source_angle = openmc.stats.Isotropic()
    source_energy = openmc.stats.Maxwell(293.6)
    source = openmc.IndependentSource(space=source_space, angle=source_angle, energy=source_energy)
    return source

@pytest.fixture
def model():
    model = openmc.Model()

    universe = create_universe()
    model.geometry = openmc.Geometry(universe)
    model.settings.run_mode = 'fixed source'
    model.settings.source = test_source()

    model.settings.batches = 10
    model.settings.particles = 100

    model.tallies = openmc.Tallies()
    tally = openmc.Tally(name='νΣf tally')
    tally.scores = ['nu-fission']
    model.tallies.append(tally)

    return model

@pytest.fixture
def model_k():
    openmc.reset_auto_ids()
    model = openmc.Model()

    universe = create_universe()
    model.geometry = openmc.Geometry(universe)
    model.settings.run_mode = 'fixed source'
    model.settings.source = test_source()
    model.settings.calculate_subcritical_k = True

    model.settings.batches = 10
    model.settings.particles = 100

    model.tallies = openmc.Tallies()
    tally = openmc.Tally(name='νΣf tally')
    tally.scores = ['nu-fission']
    model.tallies.append(tally)

    return model

def remove_xml_files():
    for f in glob.glob("*.xml"):
        os.remove(f)

def test_tally_results(model, model_k):
    """Test that fixed source run with subcritical k calculation gives the same tally results."""
    harness = PyAPITestHarness("statepoint.10.h5", model, inputs_true='inputs_true_1.dat', results_true='results_true_1.dat')
    harness.main()

    harness_k = PyAPITestHarness("statepoint.10.h5", model_k, inputs_true='inputs_true_2.dat', results_true='results_true_1.dat')
    harness_k.main()
    remove_xml_files()

def test_subcritical_k(model_k):
    """Test that subcritical k calculation gives a known result."""
    harness = PyAPITestHarness("statepoint.10.h5", model_k, inputs_true='inputs_true_3.dat', results_true='results_true_2.dat', 
                               subcritical_k_results=True)
    harness.main()
    remove_xml_files()

def test_subcritical_k_mpi_specific(model_k):
    if not config['mpi']:
        remove_xml_files()
        pytest.skip("This test is intended for MPI execution only")
        
    harness = PyAPITestHarness(
        "statepoint.10.h5", 
        model_k, 
        inputs_true='inputs_true_3.dat', 
        results_true='results_true_2.dat', 
        subcritical_k_results=True
    )
    harness.main()
    remove_xml_files()