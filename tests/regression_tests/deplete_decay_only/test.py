""" Transport-free depletion test suite """

from pathlib import Path
import shutil

import numpy as np
import pytest
import openmc
import openmc.deplete
from openmc.deplete import CoupledOperator, IndependentOperator, MicroXS


@pytest.fixture(scope="module")
def model():
    fuel = openmc.Material(name="uo2")
    fuel.add_element("U", 1, percent_type="ao", enrichment=4.25)
    fuel.add_element("O", 2)
    fuel.add_nuclide("Xe135_m1", 1)
    fuel.add_nuclide("Cs135_m1", 1)
    fuel.set_density("g/cc", 10.4)

    clad = openmc.Material(name="clad")
    clad.add_element("Zr", 1)
    clad.set_density("g/cc", 6)

    water = openmc.Material(name="water")
    water.add_element("O", 1)
    water.add_element("H", 2)
    water.set_density("g/cc", 1.0)
    water.add_s_alpha_beta("c_H_in_H2O")

    radii = [0.42, 0.45]
    fuel.volume = np.pi * radii[0] ** 2

    materials = openmc.Materials([fuel, clad, water])

    pin_surfaces = [openmc.ZCylinder(r=r) for r in radii]
    pin_univ = openmc.model.pin(pin_surfaces, materials)
    bound_box = openmc.model.RectangularPrism(1.24, 1.24, boundary_type="reflective")
    root_cell = openmc.Cell(fill=pin_univ, region=-bound_box)
    geometry = openmc.Geometry([root_cell])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 5
    settings.batches = 10

    return openmc.Model(geometry, materials, settings)

@pytest.fixture(scope="module")
def micro_xs():
    micro_xs_file = Path(__file__).parents[2] / 'micro_xs_simple.csv'
    return MicroXS.from_csv(micro_xs_file, lineterminator='\n')


@pytest.fixture(scope="module")
def chain_file():
    return Path(__file__).parents[2] / 'chain_simple_decay.xml'


@pytest.mark.parametrize("operator_type", ["coupled", "independent"])
def test_decay_only(run_in_tmpdir, operator_type, model, micro_xs, chain_file):
    """Transport free system test suite.

    """
    # Create operator
    if operator_type == "coupled":
        op = CoupledOperator(model, chain_file=chain_file)
    else:
        op = IndependentOperator(openmc.Materials([model.materials[0]]),
                                 [1e15],
                                 [micro_xs],
                                 chain_file)

    # Power and timesteps
    dt = [917.4, 2262.6]  # one Xe135_m1 half life and one Cs135_m1 half life

    # Perform simulation using the predictor algorithm
    openmc.deplete.PredictorIntegrator(op,
                                       dt,
                                       power=0.0,
                                       timestep_units='s').integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'

    # Load the reference/test results
    res_test = openmc.deplete.Results(path_test)

    _, xe135m1_atoms = res_test.get_atoms('1', 'Xe135_m1')
    _, xe135_atoms = res_test.get_atoms('1', 'Xe135')
    _, cs135m1_atoms = res_test.get_atoms('1', 'Cs135_m1')
    _, cs135_atoms = res_test.get_atoms('1', 'Cs135')

    tol = 1.0e-14
    assert xe135m1_atoms[0] == pytest.approx(xe135m1_atoms[1] * 2, rel=tol)

    # WARNING: this is generally not true as Xe135_m1 has two
    # decay modes, and Xe135 will also decay, but we've modified the depletion chain so
    # that Xe135_m1 only decays to Xe135, and that Xe135 has has no decay modes
    assert xe135_atoms[1] == pytest.approx(xe135m1_atoms[1], rel=tol)
    assert cs135m1_atoms[0] == pytest.approx(cs135m1_atoms[2] * 2, rel=tol)
    assert cs135_atoms[2] == pytest.approx(cs135m1_atoms[2], rel=tol)
