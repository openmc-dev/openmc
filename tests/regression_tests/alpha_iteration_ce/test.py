import os

import numpy as np

import openmc
from openmc.examples import pwr_pin_cell

from tests.testing_harness import PyAPITestHarness

def test_alpha_static_ce():
    model = pwr_pin_cell(alpha_mode=True)
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
