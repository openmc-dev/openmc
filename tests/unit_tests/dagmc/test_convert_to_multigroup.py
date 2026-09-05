"""Test that convert_to_multigroup works with DAGMC models without requiring
particles/batches to be set beforehand.
"""

from pathlib import Path
import pytest
import openmc
import openmc.lib

pytestmark = pytest.mark.skipif(
    not openmc.lib.feature_enabled('dagmc'),
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
        particles=10,
        overwrite_mgxs_library=True
    )

    # Verify the model was converted successfully
    assert model.settings.energy_mode == 'multi-group'


def test_convert_to_multigroup_cell_wise(run_in_tmpdir):
    """cell_wise gives each DAGMC volume its own cross sections, so two
    cells filled with the same material end up with distinct macroscopics."""
    openmc.reset_auto_ids()

    # dagmc.h5m has two fuel volumes (both "no-void fuel"), one water volume, a
    # graveyard and an implicit complement.
    u235 = openmc.Material(name="no-void fuel")
    u235.add_nuclide("U235", 1.0)
    u235.set_density("g/cm3", 11.0)
    water = openmc.Material(name="water")
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cm3", 1.0)
    water.id = 41

    dagmc_file = Path(__file__).parent / "dagmc.h5m"
    model = openmc.Model()
    model.materials = openmc.Materials([u235, water])
    model.geometry = openmc.Geometry(openmc.DAGMCUniverse(dagmc_file))
    model.settings = openmc.Settings()
    model.settings.run_mode = "fixed source"
    source = openmc.IndependentSource()
    source.energy = openmc.stats.delta_function(2.0e6)
    model.settings.source = source

    # Pre-create the library so MGXS generation/transport is skipped; this
    # exercises the per-cell material cloning for the DAGMC cells only.
    Path("mgxs.h5").touch()
    model.convert_to_multigroup(
        method="cell_wise", groups="CASMO-2", mgxs_path="mgxs.h5")

    # The three material-filled volumes (two fuel, one water) each get their own
    # cloned material with a distinct macroscopic; the void cells are skipped.
    assert model.settings.energy_mode == "multi-group"
    assert len(model.materials) == 3
    macros = [m._macroscopic for m in model.materials]
    assert len(set(macros)) == 3
