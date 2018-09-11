from tests.testing_harness import TestHarness

import pytest
pytestmark = pytest.mark.skipif(
    os.environ.get('CAD') != 'y',
    reason="CAD build is not enabled.")

def test_cad():
    harness = TestHarness('statepoint.15.h5')
    harness.main()    
