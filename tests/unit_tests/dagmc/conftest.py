from pathlib import Path

import openmc
import pytest


@pytest.fixture
def dagmc_legacy_path():
    """Path to the shared legacy DAGMC unit-test geometry."""
    return Path(__file__).parent / "dagmc.h5m"


@pytest.fixture
def dagmc_legacy_universe(dagmc_legacy_path):
    """Fresh DAGMCUniverse backed by the shared legacy test geometry."""
    return openmc.DAGMCUniverse(dagmc_legacy_path)
