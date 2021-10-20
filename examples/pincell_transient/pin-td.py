import numpy as np
import openmc
import openmc.mgxs
import openmc.kinetics as kinetics

from geometry_pin import geometry

materials_file = openmc.Materials(geometry.get_all_materials().values())

###############################################################################
# Define problem settings

settings = openmc.Settings()
settings.batches = 200
settings.inactive = 100
settings.particles = 1000
settings.output = {"tallies": False}
settings.seed = 1

# Create an initial uniform spatial source distribution over fissionable zones
lower_left = [-0.62992, -0.62992, -182.88]
upper_right = [0.62992, 0.62992, 182.88]
uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)

settings.sourcepoint = {"batches": [settings.batches], "separate": True, "write": True}

# For source convergence checks, add a mesh that can be used to calculate the
# Shannon entropy
entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = lower_left
entropy_mesh.upper_right = upper_right
entropy_mesh.dimension = [1, 1, 30]
settings.entropy_mesh = entropy_mesh

###############################################################################
# Define transient problem parameters

# Create pin cell mesh
full_pin_cell_mesh = openmc.RegularMesh()
full_pin_cell_mesh.dimension = [1, 1, 30]
full_pin_cell_mesh.lower_left = [-0.62992, -0.62992, -182.88]
full_pin_cell_mesh.width = [0.62992 * 2, 0.62992 * 2, 12.192]

# Instantiate an EnergyGroups object for the diffuion coefficients
fine_groups = openmc.mgxs.EnergyGroups()
fine_groups.group_edges = [0.0, 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.0e7]

# Instantiate an EnergyGroups object for the transient solve
energy_groups = openmc.mgxs.EnergyGroups()
energy_groups.group_edges = [0.0, 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.0e7]

# Instantiate an EnergyGroups object for one group data
one_group = openmc.mgxs.EnergyGroups()
one_group.group_edges = [fine_groups.group_edges[0], fine_groups.group_edges[-1]]

# Instantiate a clock object
t_outer = np.arange(0.0, 1.0, 5.0e-1)
clock = openmc.kinetics.Clock(dt_inner=1.0e-2, t_outer=t_outer)

# Prescribe the transient as a dictionary of densities and temperatures
transient = {}
for material in materials_file:
    transient[material.name] = {}
    for t in t_outer:
        transient[material.name][t] = {
            "density": material.density,
            "temperature": material.temperature,
        }
    if material.name == "Moderator":
        transient[material.name][0.5]["density"] = material.density * 0.9

for cell in geometry.get_all_cells().values():
    transient[cell.name] = {}
    for t in t_outer:
        transient[cell.name][t] = {"translation": None}

# Instantiate a kinetics solver object
solver = openmc.kinetics.Solver(directory="PIN_TD")
solver.num_delayed_groups = 6
solver.amplitude_mesh = full_pin_cell_mesh
solver.shape_mesh = full_pin_cell_mesh
solver.tally_mesh = full_pin_cell_mesh
solver.one_group = one_group
solver.energy_groups = energy_groups
solver.fine_groups = fine_groups
solver.tally_groups = energy_groups
solver.geometry = geometry
solver.settings = settings
solver.materials = materials_file
solver.transient = transient
solver.outer_tolerance = np.inf
solver.method = "adiabatic"
solver.multi_group = False
solver.clock = clock
solver.run_kwargs = {"threads": 1, "mpi_args": None}
solver.core_volume = np.pi * 0.39218 ** 2 * 365.76
solver.min_outer_iters = 1
solver.use_pcmfd = False
solver.use_agd = False
solver.condense_dif_coef = True
solver.chi_delayed_by_delayed_group = True
solver.chi_delayed_by_mesh = False
solver.use_pregenerated_sps = False
solver.log_file_name = "log_file.h5"

solver.solve()
