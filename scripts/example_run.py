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

# Depletion settings
settings = openmc.deplete.OpenMCSettings()
settings.power = 2.337e15*4*JOULE_PER_EV*1e6  # MeV/second cm from CASMO
settings.dt_vec = dt
settings.output_dir = 'test'

# OpenMC-delegated settings
settings.particles = 1000
settings.batches = 100
settings.inactive = 40
settings.source = openmc.Source(space=openmc.stats.Box(lower_left, upper_right))
settings.verbosity = 3

op = openmc.deplete.OpenMCOperator(geometry, settings)

# Perform simulation using the MCNPX/MCNP6 algorithm
openmc.deplete.integrator.cecm(op)
