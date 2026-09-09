import openmc
from tests.testing_harness import PyAPITestHarness


def test_weight_windows_exodus():
    # Build a minimal basic model shell to trigger initial setup parsing
    model = openmc.Model()

    # Simple fuel material
    uo2 = openmc.Material(name="uo2")
    uo2.add_nuclide("U235", 1.0)
    uo2.set_density("g/cm3", 10.0)
    model.materials.append(uo2)

    # Tiny box geometry
    box = openmc.model.RectangularPrism(
        width=0.6, height=0.9, origin=(0.3, 0.45), boundary_type="vacuum"
    )
    z_top = openmc.ZPlane(z0=0.3, boundary_type="vacuum")
    z_bot = openmc.ZPlane(z0=-0.0, boundary_type="vacuum")
    cell = openmc.Cell(fill=uo2, region=-box & +z_bot & -z_top)
    model.geometry = openmc.Geometry([cell])

    # Setup standard simulation run settings
    model.settings.batches = 3
    model.settings.inactive = 0
    model.settings.particles = 100
    model.settings.run_mode = "fixed source"
    model.settings.seed = 12345
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0.5, 0.8, 0.2))
    )
    # Configure your new Exodus feature
    model.settings.weight_windows_exodus = openmc.WeightWindowsExodus(
        file="test_out.e",  # Path to your small test mesh in this folder
        adjoint_flux_variables=["adjoint_flux_g1", "adjoint_flux_g0"],
        energy_bounds=[0.0, 0.625, 2.0e7],
    )

    # Use OpenMC's test harness to run the simulation
    harness = PyAPITestHarness("statepoint.3.h5", model)
    harness.main()
