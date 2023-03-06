""" Tests for Msr continuos class """

from pathlib import Path
import numpy as np
from math import exp
import pytest

from openmc.deplete import CoupledOperator
import openmc
from openmc.deplete import MsrContinuous
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
    geometry = openmc.Geometry([cell_f,cell_w])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)

def test_get_set(model):
    """Tests the get/set methods"""

    # create removal rates for U and Xe
    removal_rates = {'U': 0.01, 'Xe': 0.1}
    op = CoupledOperator(model, CHAIN_PATH)
    msr = MsrContinuous(op, model)

    # Test by Openmc material, material name and material id
    mat, dest_mat = [m for m in model.materials if m.depletable]
    for i in [mat, mat.name, mat.id]:
        for j in [dest_mat, dest_mat.name, dest_mat.id]:
            for key, value in removal_rates.items():
                msr.set_removal_rate(i, [key], value, destination_material=j)
                assert msr.get_removal_rate(i, key) == value
                assert msr.get_destination_material(i, key) == str(dest_mat.id)
            assert msr.get_elements(i) == removal_rates.keys()

@pytest.mark.parametrize("removal_rate_units, units_per_seconds", [
    ('1/s', 1),
    ('1/sec', 1),
    ('1/min', 60),
    ('1/minute', 60),
    ('1/h', 3600),
    ('1/hr', 3600),
    ('1/hour', 3600),
    ('1/d', 86400),
    ('1/day', 86400),
    ('1/a', 31557600.0),
    ('1/year', 31557600.0),
    ])
def test_units(removal_rate_units, units_per_seconds, model):
    """ Units testing"""
    # create removal rate Xe
    element = 'Xe'
    removal_rate = 1e-5
    op = CoupledOperator(model, CHAIN_PATH)
    msr = MsrContinuous(op, model)

    msr.set_removal_rate('f', [element], removal_rate*units_per_seconds,
                         removal_rate_units=removal_rate_units)
    assert msr.get_removal_rate('f', element) == removal_rate


def test_msr(run_in_tmpdir, model):

    """Tests msr depletion class without neither reaction rates nor decay
    but only removal rates"""

    # create removal rate for U
    element = ['U']
    removal_rate = 1e-5
    op = CoupledOperator(model, CHAIN_PATH)
    msr = MsrContinuous(op, model)
    msr.set_removal_rate('f', element, removal_rate)
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1,1], 0.0, msr_continuous = msr, timestep_units = 'd')
    integrator.integrate()

    # Get number of U238 atoms from results
    results = openmc.deplete.Results('depletion_results.h5')
    _, atoms = results.get_atoms(model.materials[0], "U238")

    # Ensure number of atoms equal removal decay
    assert atoms[1] == pytest.approx(atoms[0]*exp(-removal_rate*3600*24))
    assert atoms[2] == pytest.approx(atoms[1]*exp(-removal_rate*3600*24))
