"""Basic unit tests for openmc.deplete.Operator instantiation

Modifies and resets environment variable OPENMC_CROSS_SECTIONS
to a custom file with new depletion_chain node
"""

from os import environ
from unittest import mock
from pathlib import Path

import pytest
from openmc.deplete.abc import TransportOperator
from openmc.deplete.chain import Chain

BARE_XS_FILE = "bare_cross_sections.xml"
CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"


@pytest.fixture()
def bare_xs(run_in_tmpdir):
    """Create a very basic cross_sections file, return simple Chain.

    """

    bare_xs_contents = """<?xml version="1.0"?>
<cross_sections>
  <depletion_chain path="{}" />
</cross_sections>
""".format(CHAIN_PATH)

    with open(BARE_XS_FILE, "w") as out:
        out.write(bare_xs_contents)

    yield


class BareDepleteOperator(TransportOperator):
    """Very basic class for testing the initialization."""

    # declare abstract methods so object can be created
    def __call__(self, *args, **kwargs):
        pass

    def initial_condition(self):
        pass

    def get_results_info(self):
        pass


@mock.patch.dict(environ, {"OPENMC_CROSS_SECTIONS": BARE_XS_FILE})
def test_operator_init(bare_xs):
    """The test will set and unset environment variable OPENMC_CROSS_SECTIONS
    to point towards a temporary dummy file. This file will be removed
    at the end of the test, and only contains a
    depletion_chain node."""
    # force operator to read from OPENMC_CROSS_SECTIONS
    bare_op = BareDepleteOperator(chain_file=None)
    act_chain = bare_op.chain
    ref_chain = Chain.from_xml(CHAIN_PATH)
    assert len(act_chain) == len(ref_chain)
    for name in ref_chain.nuclide_dict:
        # compare openmc.deplete.Nuclide objects
        ref_nuc = ref_chain[name]
        act_nuc = act_chain[name]
        for prop in [
                'name', 'half_life', 'decay_energy', 'reactions',
                'decay_modes', 'yield_data', 'yield_energies',
                ]:
            assert getattr(act_nuc, prop) == getattr(ref_nuc, prop), prop
