import numpy as np
import openmc
import openmc.mgxs
import openmc.plotter
import openmc.kinetics as kinetics
from geometry_pin import geometry

materials_file = openmc.Materials(geometry.get_all_materials().values())

# Create pin cell mesh
full_pin_cell_mesh = openmc.RegularMesh()
full_pin_cell_mesh.type = 'regular'
full_pin_cell_mesh.dimension = [1,1,1]
full_pin_cell_mesh.lower_left = [-0.62992,-0.62992,-182.88]
full_pin_cell_mesh.width =[0.62992*2,0.62992*2,182.88*2]

# OpenMC simulation parameters
batches = 200
inactive = 100
particles = 1000

settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.output = {'tallies': False}

# Create an initial uniform spatial source distribution over fissionable zones
lower_left = [-0.62992, -0.62992, -182.88]
upper_right = [0.62992, 0.62992, 182.88]
uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

sourcepoint = dict()
sourcepoint['batches'] = [batches]
sourcepoint['separate'] = True
sourcepoint['write'] = True
settings_file.sourcepoint = sourcepoint

entropy_mesh = openmc.RegularMesh()
entropy_mesh.type = 'regular'
entropy_mesh.dimension = [4,4,1]
entropy_mesh.lower_left = lower_left
entropy_mesh.upper_right = upper_right
settings_file.entropy_mesh = entropy_mesh

# Instantiate an EnergyGroups object for the diffuion coefficients
fine_groups = openmc.mgxs.EnergyGroups()
fine_groups.group_edges = [0., 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.0e7]

# Instantiate an EnergyGroups object for the transient solve
energy_groups = openmc.mgxs.EnergyGroups()
energy_groups.group_edges = [0., 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.0e7]
two_groups = openmc.mgxs.EnergyGroups()
two_groups.group_edges = [0., 0.63, 1.0e7]

# Instantiate an EnergyGroups object for one group data
one_group = openmc.mgxs.EnergyGroups()
one_group.group_edges = [fine_groups.group_edges[0], fine_groups.group_edges[-1]]

# Instantiate a clock object
t_outer = np.arange(0., 1.0, 5.e-1)
clock = openmc.kinetics.Clock(start=0., end=0.5, dt_inner=1.e-2, t_outer=t_outer)


# Prescribe the transient
Transient = {}

for material in materials_file:
    MatChange = {
            material.name: {},
            }
    Transient.update(MatChange)
    for t in t_outer:
        TimeDict = {
                t: {
                    'density' : {},
                    'temperature' : {},
                    }
                }
        Transient[material.name].update(TimeDict)

for material in materials_file:
    if material.name == 'Moderator':
        Transient[material.name][0]['density'] = material.density
        Transient[material.name][0.5]['density'] = material.density*0.9
#        Transient[material.name][1.0]['density'] = material.density*0.8
#        Transient[material.name][1.5]['density'] = material.density*0.9
#        Transient[material.name][2.0]['density'] = material.density
        for t in t_outer:
            Transient[material.name][t]['temperature'] = material.temperature
    else:
        for t in t_outer:
            Transient[material.name][t]['density'] = material.density
            Transient[material.name][t]['temperature'] = material.temperature

#Instantiate a kinetics solver object
solver = openmc.kinetics.Solver(directory='PIN_TD_ADIABATIC')
solver.num_delayed_groups           = 6
solver.amplitude_mesh               = full_pin_cell_mesh
solver.shape_mesh                   = full_pin_cell_mesh
solver.tally_mesh                   = full_pin_cell_mesh
solver.one_group                    = one_group
solver.energy_groups                = energy_groups
solver.fine_groups                  = fine_groups
solver.tally_groups                 = energy_groups
solver.geometry                     = geometry
solver.settings_file                = settings_file
solver.materials_file               = materials_file
solver.transient                    = Transient
solver.inner_tolerance              = np.inf
solver.outer_tolerance              = np.inf
solver.method                       = 'ADIABATIC'
solver.multi_group                  = False
solver.clock                        = clock
solver.mpi_procs                    = 1
solver.threads                      = 18
solver.core_volume                  = np.pi*0.39218**2*365.76
solver.constant_seed                = True
solver.seed                         = 1
solver.min_outer_iters              = 1
solver.use_pcmfd                    = False
solver.use_agd                      = False
solver.condense_dif_coef            = True
solver.chi_delayed_by_delayed_group = True
solver.chi_delayed_by_mesh          = False
solver.use_pregenerated_sps         = False
solver.run_on_cluster               = False
solver.log_file_name                = 'log_file.h5'

solver.solve()
