import openmc
import pytest


def test_source_rejection_fraction(run_in_tmpdir):
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    w = 0.25
    rpp1 = openmc.model.RectangularParallelepiped(-w/2, w/2, -w/2, w/2, -w/2, w/2)
    rpp2 = openmc.model.RectangularParallelepiped(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=mat, region=-rpp1)
    cell2 = openmc.Cell(region=+rpp1 & -rpp2)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2])

    # Create a box source over a 1 cm³ volume that is constrained to the source
    # cell of volume (0.25 cm)³ = 0.0125 cm³, which means the default rejection
    # fraction of 0.05 won't work
    model.settings.particles = 100000
    model.settings.batches = 1
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(*(-rpp2).bounding_box),
        constraints={'domains': [cell1]}
    )
    with pytest.raises(RuntimeError):
        model.run()

    # With a source rejection fraction below 0.0125, the simulation should run
    model.settings.source_rejection_fraction = 0.005
    model.run()
