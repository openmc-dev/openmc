"""An example file showing how to run a simulation."""

import numpy as np
import openmc
from openmc.data import JOULE_PER_EV
import openmc.deplete

import example_geometry

# Load geometry from example
geometry, lower_left, upper_right = example_geometry.generate_problem()

# Create dt vector for 5.5 months with 15 day timesteps
dt1 = 15*24*60*60  # 15 days
dt2 = 5.5*30*24*60*60  # 5.5 months
N = np.floor(dt2/dt1)
dt = np.repeat([dt1], N)

# Power for simulation
power = 2.337e15*4*JOULE_PER_EV*1e6  # MeV/second cm from CASMO

# OpenMC settings
settings = openmc.Settings()
settings.particles = 1000
settings.batches = 100
settings.inactive = 40
settings.source = openmc.Source(space=openmc.stats.Box(lower_left, upper_right))

op = openmc.deplete.Operator(geometry, settings)
op.output_dir = 'test'

# Perform simulation using the MCNPX/MCNP6 algorithm
openmc.deplete.integrator.cecm(op, dt, power)
