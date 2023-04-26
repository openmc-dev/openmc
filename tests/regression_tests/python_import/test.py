import sys
import openmc


def test_matplotlib_presence():
    """Checks that remains a deferred import"""
    assert 'matplotlib.pyplot' not in sys.modules
