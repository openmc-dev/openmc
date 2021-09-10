import openmc
import math

###############################################################################
#                             Design Parameters
###############################################################################

thermal_cutoff = 4.0
mpi_procs = 4

###############################################################################
#                          Create settings file
###############################################################################

# Instantiate a Settings object
settings_file               = openmc.Settings()
settings_file.batches       = 200
settings_file.inactive      = 100
settings_file.particles     = 1000
settings_file.output        = {'tallies': True}

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-32.13, -10.71, -64.26,
           10.71,  32.13,  64.26]

uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

settings_file.entropy_lower_left  = bounds[:3]
settings_file.entropy_upper_right = bounds[3:]
settings_file.entropy_dimension   = [34,34,1]
