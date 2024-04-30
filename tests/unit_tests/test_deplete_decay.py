from pathlib import Path

import openmc.deplete
import numpy as np
import pytest


def test_deplete_decay_products(run_in_tmpdir):
    # Create chain file with H1, He4, and Li5
    with open('test_chain.xml', 'w') as chain_file:
        chain_file.write("""
<depletion_chain>
  <nuclide name="H1" reactions="0">
  </nuclide>
  <nuclide name="He4" reactions="0"/>
  <nuclide name="Li5" half_life="3.06868e-22" decay_modes="1" decay_energy="1965000.0" reactions="0">
    <decay type="alpha" target="H1" branching_ratio="1.0"/>
  </nuclide>
</depletion_chain>
        """)

    # Create MicroXS object with no cross sections
    micro_xs = openmc.deplete.MicroXS(np.empty((0, 0)), [], [])

    # Create depletion operator with no reactions
    op = openmc.deplete.IndependentOperator.from_nuclides(
        volume=1.0,
        nuclides={'Li5': 1.0},
        flux=0.0,
        micro_xs=micro_xs,
        chain_file='test_chain.xml',
        normalization_mode='source-rate'
    )

    # Create time-integrator and integrate
    integrator = openmc.deplete.PredictorIntegrator(
        op, timesteps=[1.0], source_rates=[0.0], timestep_units='d'
    )
    integrator.integrate(final_step=False)

    # Get concentration of H1 and He4
    results = openmc.deplete.Results('depletion_results.h5')
    _, h1 = results.get_atoms("1", "H1")
    _, he4 = results.get_atoms("1", "He4")

    # Since we started with 1e24 atoms of Li5, we should have 1e24 atoms of both
    # H1 and He4
    assert h1[1] == pytest.approx(1e24)
    assert he4[1] == pytest.approx(1e24)


def test_deplete_decay_step_fissionable(run_in_tmpdir):
    """Ensures that power is not computed in zero power cases with
    fissionable material present. This tests decay calculations without
    power, although this specific example does not exhibit any decay.

    Proves github issue #2963 is fixed
    """

    # Set up a pure decay operator
    micro_xs = openmc.deplete.MicroXS(np.empty((0, 0)), [], [])
    mat = openmc.Material()
    mat.name = 'I do not decay.'
    mat.add_nuclide('U238', 1.0, 'ao')
    mat.volume = 10.0
    mat.set_density('g/cc', 1.0)
    original_atoms = mat.get_nuclide_atoms()['U238']

    mats = openmc.Materials([mat])
    op = openmc.deplete.IndependentOperator(
        mats, [1.0], [micro_xs], Path(__file__).parents[1] / "chain_simple.xml")

    # Create time integrator and integrate
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1.0], power=[0.0], timestep_units='s'
    )
    integrator.integrate()

    # Get concentration of U238. It should be unchanged since this chain has no U238 decay.
    results = openmc.deplete.Results('depletion_results.h5')
    _, u238 = results.get_atoms("1", "U238")

    assert u238[1] == pytest.approx(original_atoms)
