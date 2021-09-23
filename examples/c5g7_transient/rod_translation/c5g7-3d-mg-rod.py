import numpy as np
import os
import openmc
import openmc.mgxs
import openmc.kinetics as kinetics

from geometry_mg import materials, surfaces, universes, cells, lattices, geometry, mgxs_lib_file
from mgxs_lib import mgxs_data
materials_file = openmc.Materials(geometry.get_all_materials().values())
os.environ['OPENMC_MG_CROSS_SECTIONS'] = 'mgxs.h5'

###############################################################################
# Define problem settings

settings = openmc.Settings()
settings.batches = 110
settings.inactive = 40
settings.particles = 1000
settings.output = {'tallies': False}
settings.seed = 1
settings.energy_mode = 'multi-group'

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-32.13, -10.71, -64.26, 10.71,  32.13,  64.26]
entropy_bounds = [-32.13, -10.71, -85.68, 10.71,  32.13,  85.68]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)

settings.sourcepoint = {
    'batches': [settings.batches],
    'separate': True,
    'write': True
}

# For source convergence checks, add a mesh that can be used to calculate the
# Shannon entropy
entropy_mesh = openmc.RegularMesh()
entropy_mesh.dimension = [4,4,32]
entropy_mesh.lower_left  = entropy_bounds[:3]
entropy_mesh.upper_right = entropy_bounds[3:]
settings.entropy_mesh = entropy_mesh

###############################################################################
# Define transient problem parameters

full_pin_cell_mesh = openmc.RegularMesh()
full_pin_cell_mesh.dimension = [51,51,3]
full_pin_cell_mesh.lower_left  = [-32.13, -32.13, -64.26]
full_pin_cell_mesh.upper_right = [ 32.13,  32.13,  64.26]

full_assembly_mesh = openmc.RegularMesh()
full_assembly_mesh.dimension = [3,3,3]
full_assembly_mesh.lower_left  = [-32.13, -32.13, -64.26]
full_assembly_mesh.upper_right = [ 32.13,  32.13,  64.26]

# Instantiate a 50-group EnergyGroups object
fine_groups = openmc.mgxs.EnergyGroups()
fine_groups.group_edges = [0., 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.0e7]

# Instantiate a 2-group EnergyGroups object
energy_groups = openmc.mgxs.EnergyGroups()
energy_groups.group_edges = [0., 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.0e7]

# Instantiate a 1-group EnergyGroups object
one_group = openmc.mgxs.EnergyGroups()
one_group.group_edges = [fine_groups.group_edges[0], fine_groups.group_edges[-1]]

# Instantiate a clock object
t_outer = np.arange(0., 1.5, 5.e-1)
clock = openmc.kinetics.Clock(dt_inner=1.e-2, t_outer=t_outer)

# Prescribe the transient as a dictionary of densities and temperatures 
transient = {}
for material in materials_file:
    transient[material.name] = {}
    for t in t_outer:
        transient[material.name][t] = {
            'density': material.density,
            'temperature': material.temperature
        }
for cell in geometry.get_all_cells().values():
	transient[cell.name] = {}
	for t in t_outer:
		transient[cell.name][t] = {
		     'translation': None
		}
transient['Control Rod Base Bank 1'][0.0]['translation'] = [0., 0., 64.26]
transient['Control Rod Base Bank 1'][0.5]['translation'] = [0., 0., 21.42]
transient['Control Rod Base Bank 1'][1.0]['translation'] = [0., 0., 64.26]

# Instantiate a kinetics solver object
solver = openmc.kinetics.Solver(directory='C5G7_3D_MG')
solver.num_delayed_groups           = 6
solver.amplitude_mesh               = full_assembly_mesh
solver.shape_mesh                   = full_pin_cell_mesh
solver.tally_mesh                   = full_pin_cell_mesh
solver.one_group                    = one_group
solver.energy_groups                = energy_groups
solver.fine_groups                  = energy_groups
solver.tally_groups                 = energy_groups
solver.geometry                     = geometry
solver.settings                     = settings
solver.materials                    = materials_file
solver.transient                    = transient
solver.outer_tolerance              = np.inf
solver.mgxs_lib                     = mgxs_lib_file
solver.method                       = 'adiabatic'
solver.multi_group                  = True
solver.clock                        = clock
solver.run_kwargs                   = {'threads': 1, 'mpi_args': None}
solver.core_volume                  = 42.84 * 42.84 * 128.52
solver.min_outer_iters              = 1
solver.use_pcmfd                    = False
solver.use_agd                      = False
solver.condense_dif_coef            = True
solver.chi_delayed_by_delayed_group = True
solver.chi_delayed_by_mesh          = False
solver.use_pregenerated_sps         = False
solver.log_file_name                = 'log_file.h5'

# Run OpenMC
solver.solve()