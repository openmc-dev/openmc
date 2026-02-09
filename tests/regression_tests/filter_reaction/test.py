from openmc.examples import pwr_pin_cell
from tests.testing_harness import PyAPITestHarness

import openmc


def test_filter_reaction():
    model = pwr_pin_cell()

    # Create a tally with reaction filter
    tally = openmc.Tally()
    tally.filters = [openmc.ReactionFilter(
        ['(n,elastic)', '(n,fission)', '(n,gamma)']
    )]
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    # Reduce particles for faster testing
    model.settings.particles = 1000
    model.settings.batches = 10

    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
