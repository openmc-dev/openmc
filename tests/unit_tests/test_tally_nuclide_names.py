"""A tally nuclide bin names a nuclide.

Anything else has to be reported rather than terminating on an uncaught
exception, which is what a name the data library cannot resolve used to do.
"""

import openmc
import pytest


def model_with_tally_nuclides(nuclides):
    openmc.reset_auto_ids()

    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 10.0)

    sphere = openmc.Sphere(r=5.0, boundary_type='vacuum')
    model = openmc.Model()
    model.geometry = openmc.Geometry([openmc.Cell(fill=mat, region=-sphere)])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 100
    model.settings.batches = 1
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())

    tally = openmc.Tally()
    tally.nuclides = nuclides
    tally.scores = ['total']
    model.tallies = openmc.Tallies([tally])
    return model


@pytest.mark.parametrize('name', ['Al', 'U', 'all', 'not-a-nuclide'])
def test_unresolvable_nuclide_is_reported(run_in_tmpdir, name):
    """An element symbol is not a nuclide, and neither is a typo."""
    model = model_with_tally_nuclides([name])
    with pytest.raises(RuntimeError, match='Could not add nuclide'):
        model.run()


def test_element_symbol_of_a_present_element(run_in_tmpdir):
    """Reported even when the element itself is in the material.

    The bin is matched against nuclide names, so 'U' does not resolve in a
    material made of U235.
    """
    model = model_with_tally_nuclides(['U'])
    with pytest.raises(RuntimeError, match='Could not add nuclide'):
        model.run()


def test_total_and_nuclide_bins_together(run_in_tmpdir):
    """The total bin stays valid alongside real nuclides."""
    model = model_with_tally_nuclides(['total', 'U235'])
    model.run()
