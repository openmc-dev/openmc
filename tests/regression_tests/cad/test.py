from tests.testing_harness import CADTestHarness
import os
import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get('CAD') != 'y',
    reason="CAD build is not enabled.")

def test_cad():
    harness = CADTestHarness('statepoint.15.h5')
    harness.main()    
