import pytest
from sandbox.tokomak_analysis import Config,build_model


@pytest.fixture(scope='module')
def assembly():return build_model(Config(),tallies=False)
