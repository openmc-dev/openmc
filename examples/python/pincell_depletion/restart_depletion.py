import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 100
inactive = 10
particles = 1000

# Depletion simulation parameters
time_step = 1*24*60*60 # s
final_time = 5*24*60*60 # s
time_steps = np.full(final_time // time_step, time_step)

chain_file = './chain_simple.xml'
power = 174 # W/cm, for 2D simulations only (use W for 3D)

###############################################################################
#                      Load previous simulation results
###############################################################################

# Load geometry from statepoint
statepoint = 'statepoint.100.h5'
sp = openmc.StatePoint(statepoint)
geometry = sp.summary.geometry

# Close statepoint and summary files to be able to write over them
sp.close()

# Load previous depletion results
previous_results = openmc.deplete.ResultsList("depletion_results.h5")

###############################################################################
#                      Transport calculation settings
###############################################################################

# Instantiate a Settings object, set all runtime parameters
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

entropy_mesh = openmc.Mesh()
entropy_mesh.lower_left = [-0.39218, -0.39218, -1.e50]
entropy_mesh.upper_right = [0.39218, 0.39218, 1.e50]
entropy_mesh.dimension = [10, 10, 1]
settings_file.entropy_mesh = entropy_mesh

###############################################################################
#                   Initialize and run depletion calculation
###############################################################################

op = openmc.deplete.Operator(geometry, settings_file, chain_file, \
                                                  previous_results)

# Perform simulation using the predictor algorithm
openmc.deplete.integrator.predictor(op, time_steps, power)

###############################################################################
#                    Read depletion calculation results
###############################################################################

# Open results file
results = openmc.deplete.ResultsList("depletion_results.h5")

# Obtain K_eff as a function of time
time, keff = results.get_eigenvalue()
                                                                                                                                                                                                                                                                                                                                                                                                                                                        
# Plot eigenvalue as a function of time
plt.figure()
plt.plot(time/24/60/60, keff, label="K-effective")
plt.xlabel("Time (days)")
plt.ylabel("Keff")
plt.show()
