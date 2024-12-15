"""Basic unit tests for openmc.deplete.Operator instantiation

"""

from pathlib import Path

from openmc.deplete.abc import TransportOperator
from openmc.deplete.chain import Chain

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"


class BareDepleteOperator(TransportOperator):
    """Very basic class for testing the initialization."""

    @staticmethod
    def __call__(*args, **kwargs):
        pass

    @staticmethod
    def initial_condition():
        pass

    @staticmethod
    def get_results_info():
        pass

    @staticmethod
    def write_bos_data():
        pass


def test_operator_init():
    """The test uses a temporary dummy chain. This file will be removed
    at the end of the test, and only contains a depletion_chain node."""
    bare_op = BareDepleteOperator(CHAIN_PATH)
    act_chain = bare_op.chain
    ref_chain = Chain.from_xml(CHAIN_PATH)
    assert len(act_chain) == len(ref_chain)
    for name in ref_chain.nuclide_dict:
        # compare openmc.deplete.Nuclide objects
        ref_nuc = ref_chain[name]
        act_nuc = act_chain[name]
        for prop in [
            "name",
            "half_life",
            "decay_energy",
            "reactions",
            "decay_modes",
            "yield_data",
            "yield_energies",
        ]:
            assert getattr(act_nuc, prop) == getattr(ref_nuc, prop), prop


def test_operator_fiss_q():
    """Make sure fission q values can be set"""
    new_q = {"U235": 2.0e8, "U238": 2.0e8, "U234": 5.0e7}
    operator = BareDepleteOperator(chain_file=CHAIN_PATH, fission_q=new_q)
    mod_chain = operator.chain
    for name, q in new_q.items():
        chain_nuc = mod_chain[name]
        for rx in chain_nuc.reactions:
            if rx.type == "fission":
                assert rx.Q == q
