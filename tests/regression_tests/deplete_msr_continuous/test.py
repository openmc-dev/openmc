""" MSR depletion test suite """

from pathlib import Path

import numpy as np
import pytest
import openmc
import openmc.deplete
from openmc.deplete import CoupledOperator
from openmc.deplete import MsrContinuous

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
    surf_w = openmc.Sphere(r=radii[1], boundary_type='reflective')
    cell_f = openmc.Cell(fill=f, region=-surf_f)
    cell_w = openmc.Cell(fill=w, region=+surf_f & -surf_w)
    geometry = openmc.Geometry([cell_f,cell_w])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)

@pytest.mark.parametrize("removal, feed, power, ref_result", [
    (True, False, 0.0, 'no_depletion_only_removal'),
    (False, True, 0.0, 'no_depletion_only_feed'),
    (True, False, 174.0, 'depletion_with_removal'),
    (False, True, 174.0, 'depletion_with_feed'),
    (True, True, 0.0, 'no_depletion_with_feed_and_removal'),
    (True, True, 174.0, 'depletion_with_feed_and_removal'),
    ])
def test_msr(run_in_tmpdir, model, removal, feed, power, ref_result):
    """Tests msr depletion class with removal rates"""

    chain_file = Path(__file__).parents[2] / 'chain_simple.xml'

    # Defining removing and feeding elements, using the same removal rate
    removing_element = ['Xe']
    feeding_element = ['U']
    removal_rate = 1e-5

    op = CoupledOperator(model, chain_file)
    msr = MsrContinuous(op, model)
    if removal:
        msr.set_removal_rate('f', removing_element, removal_rate,
                            destination_material=None)
    if feed:
        msr.set_removal_rate('f', feeding_element, removal_rate,
                            destination_material='w')
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1], power, msr_continuous = msr, timestep_units = 'd')
    integrator.integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'
    path_reference = Path(__file__).with_name(f'ref_msr_{ref_result}.h5')

    # Load the reference/test results
    res_ref = openmc.deplete.Results(path_reference)
    res_test = openmc.deplete.Results(path_test)

    for mat in res_test[0].rates[0].index_mat:
        for nuc in res_test[0].rates[0].index_nuc:
            for rx in res_test[0].rates[0].index_rx:
                y_test = res_test.get_reaction_rate(mat, nuc, rx)[1] / \
                    res_test.get_atoms(mat, nuc)[1]
                y_ref = res_ref.get_reaction_rate(mat, nuc, rx)[1] / \
                    res_ref.get_atoms(mat, nuc)[1]
                assert y_test == pytest.approx(y_ref, abs=1e-6)
