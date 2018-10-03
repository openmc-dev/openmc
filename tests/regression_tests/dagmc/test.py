from tests.testing_harness import TestHarness
import os
import pytest
import openmc

pytestmark = pytest.mark.skipif(
    not openmc.capi.dagmc_enabled,
    reason="DAGMC CAD geometry is not enabled.")

def test_dagmc():
    harness = TestHarness('statepoint.5.h5')
    harness.main()    
