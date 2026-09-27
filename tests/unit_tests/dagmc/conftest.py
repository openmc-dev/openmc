from pathlib import Path

import pytest


@pytest.fixture
def dagmc_legacy_path() -> Path:
    """Return the shared legacy DAGMC model used by unit tests."""
    path = Path(__file__).parent / "dagmc.h5m"
    assert path.exists()
    return path
