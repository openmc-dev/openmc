from tests.testing_harness import CADTestHarness
import os
import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get('DAGMC') != 'y',
    reason="DAGMC CAD geometry is not enabled.")

def test_cad():
    harness = CADTestHarness('statepoint.5.h5')
    harness.main()    
