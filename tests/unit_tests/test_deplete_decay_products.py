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
