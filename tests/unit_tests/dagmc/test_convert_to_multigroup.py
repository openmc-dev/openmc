"""Test that convert_to_multigroup works with DAGMC models without requiring
particles/batches to be set beforehand.
"""

from pathlib import Path
import pytest
import openmc
import openmc.lib

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")


def test_convert_to_multigroup_without_particles_batches(run_in_tmpdir):
    """Test that convert_to_multigroup works with DAGMC model without
    setting particles/batches beforehand."""
    openmc.reset_auto_ids()
    
    mat = openmc.Material(name="mat")
    mat.add_nuclide("Fe56", 1.0)
    mat.set_density("g/cm3", 7.0)

    # Use minimal tetrahedral DAGMC file
    dagmc_file = Path(__file__).parent / "dagmc_tetrahedral_no_graveyard.h5m"
    dagmc_univ = openmc.DAGMCUniverse(dagmc_file, auto_geom_ids=True)
    bound_dagmc_univ = dagmc_univ.bounded_universe(padding_distance=1)

    # Create model WITHOUT setting particles or batches
    model = openmc.Model()
    model.materials = openmc.Materials([mat])
    model.geometry = openmc.Geometry(bound_dagmc_univ)
    model.settings = openmc.Settings()  # Note: no particles or batches set!

    model.settings.run_mode = 'fixed source'
    
    # Create a point source
    my_source = openmc.IndependentSource()
    my_source.space = openmc.stats.Point((0.25, 0.25, 0.25))
    my_source.energy = openmc.stats.delta_function(14e6)
    model.settings.source = my_source

    # This should work without requiring particles/batches to be set
    # convert_to_multigroup handles initialization internally using non-transport mode
    model.convert_to_multigroup(
        method='material_wise',
        groups='CASMO-2',
        nparticles=10,
        overwrite_mgxs_library=True
    )

    # Verify the model was converted successfully
    assert model.settings.energy_mode == 'multi-group'
