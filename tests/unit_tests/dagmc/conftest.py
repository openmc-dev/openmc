
import pytest

@pytest.fixture(scope='module', autouse=True)
def setup_dagmc_unit_test(request):

    # Change to test directory
    olddir = request.fspath.dirpath().chdir()
    try:
        yield
    finally:
        olddir.chdir()
