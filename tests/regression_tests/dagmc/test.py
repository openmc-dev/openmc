from tests.testing_harness import TestHarness
import os
import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get('DAGMC') != 'y',
    reason="DAGMC CAD geometry is not enabled.")

def test_dagmc():
    harness = TestHarness('statepoint.5.h5')
    harness.main()    
