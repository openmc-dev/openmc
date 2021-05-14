import numpy as np
import os
import openmc
import openmc.mgxs
import openmc.kinetics as kinetics

from geometry_2d_mg import materials, surfaces, universes, cells, lattices, geometry, mgxs_lib_file
from mgxs_lib import mgxs_data
materials_file = openmc.Materials(geometry.get_all_materials().values())
os.environ['OPENMC_MG_CROSS_SECTIONS'] = 'mgxs.h5'

###############################################################################
# Define problem settings

settings_file = openmc.Settings()
settings_file.batches = 100
settings_file.inactive = 40
settings_file.particles = 1000
settings_file.output = {'tallies': False}
settings_file.seed = 1
settings_file.energy_mode = 'multi-group'

# Create an initial uniform spatial source distribution over fissionable zones
source_bounds  = [-32.13, -10.71, -64.26, 10.71,  32.13,  64.26]
uniform_dist = openmc.stats.Box(source_bounds[:3], source_bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

sourcepoint = dict()
sourcepoint['batches'] = [settings_file.batches]
sourcepoint['separate'] = True
sourcepoint['write'] = True
settings_file.sourcepoint = sourcepoint

# For source convergence checks, add a mesh that can be used to calculate the
# Shannon entropy
entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left  = source_bounds[:3]
entropy_mesh.upper_right = source_bounds[3:]
entropy_mesh.dimension = [4,4,1]
settings_file.entropy_mesh = entropy_mesh

###############################################################################
# Define transient problem parameters

# Create pin cell mesh
full_pin_cell_mesh = openmc.RegularMesh()
full_pin_cell_mesh.dimension = [51,51,1]
full_pin_cell_mesh.lower_left  = [-32.13, -32.13, -64.26]
full_pin_cell_mesh.width = [1.26,  1.26,  128.52]

# Create assembly mesh
full_assembly_mesh = openmc.RegularMesh()
full_assembly_mesh.dimension = [3,3,1]
full_assembly_mesh.lower_left  = [-32.13, -32.13, -64.26]
full_assembly_mesh.width = [ 21.42,  21.42,  128.52]

# Instantiate an EnergyGroups object for the diffuion coefficients
fine_groups = openmc.mgxs.EnergyGroups()
fine_groups.group_edges = [0., 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.0e7]

# Instantiate an EnergyGroups object for the transient solve
energy_groups = openmc.mgxs.EnergyGroups()
energy_groups.group_edges = [0., 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.0e7]

# Instantiate an EnergyGroups object for one group data
one_group = openmc.mgxs.EnergyGroups()
one_group.group_edges = [fine_groups.group_edges[0], fine_groups.group_edges[-1]]

# Instantiate a clock object
t_outer = np.arange(0., 1.0, 5.e-1)
clock = openmc.kinetics.Clock(dt_inner=1.e-2, t_outer=t_outer)

# Prescribe the transient as a dictionary of densities and temperatures 
transient = {}

# Create entries for each material 
for material in materials_file:
    MatChange = {
        material.name: {},
        }

    transient.update(MatChange)
    for t in t_outer:
        time_dict = {
            t: {
                'density' : {},
                'temperature' : {},
                }
            }
        transient[material.name].update(time_dict)

# Fill the entries with the desired values
for material in materials_file:
    for t in t_outer:
        transient[material.name][t]['density'] = material.density
        transient[material.name][t]['temperature'] = material.temperature

for bank in range(1,2):
    name = 'Moderator Bank {}'.format(bank)
    transient[name][0]['density'] = materials[name].density
    transient[name][0.5]['density'] = materials[name].density*0.9
    for t in t_outer:
        transient[name][t]['temperature'] = materials[name].temperature

# Instantiate a kinetics solver object
solver = openmc.kinetics.Solver(directory='C5G7_TD_MG')
solver.num_delayed_groups           = 6
solver.amplitude_mesh               = full_assembly_mesh
solver.shape_mesh                   = full_pin_cell_mesh
solver.tally_mesh                   = full_pin_cell_mesh
solver.one_group                    = one_group
solver.energy_groups                = energy_groups
solver.fine_groups                  = fine_groups
solver.tally_groups                 = energy_groups
solver.geometry                     = geometry
solver.settings                     = settings_file
solver.materials                    = materials_file
solver.transient                    = transient
solver.inner_tolerance              = np.inf
solver.outer_tolerance              = np.inf
solver.mgxs_lib                     = mgxs_lib_file
solver.method                       = 'ADIABATIC'
solver.multi_group                  = True
solver.clock                        = clock
solver.run_kwargs                   = {'threads':1, 'mpi_args':None}
solver.core_volume                  = 42.84 * 42.84 * 128.52
solver.min_outer_iters              = 1
solver.use_pcmfd                    = False
solver.use_agd                      = False
solver.condense_dif_coef            = True
solver.chi_delayed_by_delayed_group = True
solver.chi_delayed_by_mesh          = False
solver.use_pregenerated_sps         = False
solver.job_file                     = 'job_broadwell.pbs'
solver.log_file_name                = 'log_file.h5'

# Solve transient problem
solver.solve()
