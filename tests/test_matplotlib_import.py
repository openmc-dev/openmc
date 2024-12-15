import sys
import openmc


def test_matplotlib_presence():
    """Checks that matplotlib remains a deferred import"""
    assert "matplotlib" not in sys.modules
