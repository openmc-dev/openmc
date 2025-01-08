import sys
import openmc


def test_matplotlib_presence():
    """Checks that matplotlib remains a deferred import"""
    assert 'matplotlib' not in sys.modules


def test_openmc_lib_presence():
    """Checks that openmc.lib remains a deferred import"""
    assert 'openmc.lib' not in sys.modules
