import openmc

from tests.testing_harness import PyAPITestHarness


def test_filter_reaction():
    model = openmc.Model()

    m = openmc.Material()
    m.set_density('g/cm3', 10.0)
    m.add_nuclide('U235', 1.0)
    model.materials.append(m)

    s = openmc.Sphere(r=100.0, boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=-s)
    model.geometry = openmc.Geometry([c])

    # Create a tally with reaction filter
    tally = openmc.Tally()
    tally.filters = [openmc.ReactionFilter(
        ['(n,elastic)', '(n,2n)', '(n,fission)', '(n,gamma)', 'total']
    )]
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    # Reduce particles for faster testing
    model.settings.particles = 1000
    model.settings.batches = 5

    harness = PyAPITestHarness('statepoint.5.h5', model)
    harness.main()
