""" Tests for TransferRates class """

from pathlib import Path
from math import exp

import pytest
import numpy as np

import openmc
from openmc.deplete import CoupledOperator
from openmc.deplete.transfer_rates import _TransferRates
from openmc.deplete.abc import (_SECONDS_PER_MINUTE, _SECONDS_PER_HOUR,
                                _SECONDS_PER_DAY, _SECONDS_PER_JULIAN_YEAR)

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"

@pytest.fixture
def model():
    f = openmc.Material(name="f")
    f.add_element("U", 1, percent_type="ao", enrichment=4.25)
    f.add_element("O", 2)
    f.set_density("g/cc", 10.4)

    w = openmc.Material(name="w")
    w.add_element("O", 1)
    w.add_element("H", 2)
    w.set_density("g/cc", 1.0)
    w.depletable = True

    radii = [0.42, 0.45]
    f.volume = np.pi * radii[0] ** 2
    w.volume = np.pi * (radii[1]**2 - radii[0]**2)
    materials = openmc.Materials([f, w])

    surf_f = openmc.Sphere(r=radii[0])
    surf_w = openmc.Sphere(r=radii[1], boundary_type='vacuum')
    cell_f = openmc.Cell(fill=f, region=-surf_f)
    cell_w = openmc.Cell(fill=w, region=+surf_f & -surf_w)
    geometry = openmc.Geometry([cell_f, cell_w])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)


def test_get_set(model):
    """Tests the get/set methods"""

    # create transfer rates for U and Xe
    transfer_rates = {'U': 0.01, 'Xe': 0.1}
    op = CoupledOperator(model, CHAIN_PATH)
    transfer = _TransferRates(op, model)

    # Test by Openmc material, material name and material id
    material, dest_material = [m for m in model.materials if m.depletable]

    for material_input in [material, material.name, material.id]:
        for dest_material_input in [dest_material, dest_material.name,
                                    dest_material.id]:
            for element, transfer_rate in transfer_rates.items():
                transfer.set_transfer_rate(material_input, [element], transfer_rate,
                                 destination_material=dest_material_input)
                assert transfer.get_transfer_rate(material_input,
                                            element) == transfer_rate
                assert transfer.get_destination_material(material_input,
                                                    element) == str(dest_material.id)
            assert transfer.get_elements(material_input) == transfer_rates.keys()


@pytest.mark.parametrize("transfer_rate_units, unit_conv", [
    ('1/s', 1),
    ('1/sec', 1),
    ('1/min', _SECONDS_PER_MINUTE),
    ('1/minute', _SECONDS_PER_MINUTE),
    ('1/h', _SECONDS_PER_HOUR),
    ('1/hr', _SECONDS_PER_HOUR),
    ('1/hour', _SECONDS_PER_HOUR),
    ('1/d', _SECONDS_PER_DAY),
    ('1/day', _SECONDS_PER_DAY),
    ('1/a', _SECONDS_PER_JULIAN_YEAR),
    ('1/year', _SECONDS_PER_JULIAN_YEAR),
    ])
def test_units(transfer_rate_units, unit_conv, model):
    """ Units testing"""
    # create transfer rate Xe
    element = 'Xe'
    transfer_rate = 1e-5
    op = CoupledOperator(model, CHAIN_PATH)
    transfer = _TransferRates(op, model)

    transfer.set_transfer_rate('f', [element], transfer_rate * unit_conv,
                         transfer_rate_units=transfer_rate_units)
    assert transfer.get_transfer_rate('f', element) == transfer_rate


def test_transfer(run_in_tmpdir, model):
    """Tests transfer depletion class without neither reaction rates nor decay
    but only transfer rates"""

    # create transfer rate for U
    element = ['U']
    transfer_rate = 1e-5
    op = CoupledOperator(model, CHAIN_PATH)
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1,1], 0.0, timestep_units = 'd')
    integrator.add_transfer_rate('f', element, transfer_rate)
    integrator.integrate()

    # Get number of U238 atoms from results
    results = openmc.deplete.Results('depletion_results.h5')
    _, atoms = results.get_atoms(model.materials[0], "U238")

    # Ensure number of atoms equal transfer decay
    assert atoms[1] == pytest.approx(atoms[0]*exp(-transfer_rate*3600*24))
    assert atoms[2] == pytest.approx(atoms[1]*exp(-transfer_rate*3600*24))
