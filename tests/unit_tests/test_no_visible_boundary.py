import openmc


def test_no_visible_boundary(run_in_tmpdir):
    copper = openmc.Material()
    copper.add_nuclide("Cu63", 1.0)
    copper.set_density("g/cm3", 0.3)
    air = openmc.Material()
    air.add_nuclide("N14", 1.0)
    air.set_density("g/cm3", 0.0012)

    # Create a simple model of a neutron source directly impinging on a thin
    # disc of copper. Neutrons leaving the back of the disc see no surfaces in
    # front of them.
    disc = openmc.model.RightCircularCylinder((0.0, 0.0, 1.0), 0.1, 1.2)
    box = openmc.model.RectangularPrism(width=10, height=10, boundary_type="vacuum")
    c1 = openmc.Cell(fill=copper, region=-disc)
    c2 = openmc.Cell(fill=air, region=+disc & -box)
    model = openmc.Model()
    model.geometry = openmc.Geometry([c1, c2])
    model.settings.run_mode = "fixed source"
    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(), angle=openmc.stats.Monodirectional((0.0, 0.0, 1.0))
    )

    # Run model to ensure it doesn't segfault
    model.run()
