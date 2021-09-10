import openmc
import numpy as np
from geometry import materials, surfaces, universes, cells, lattices, geometry

###############################################################################
#                              general tallies
###############################################################################

tallies = {}

# Create energy group structures
seven_groups = openmc.EnergyFilter([0., 0.13, 0.63, 4.1, 55.6, 9.2e3, 1.36e6, 1.e7])
two_groups = openmc.EnergyFilter([0., 4.1, 1.e7])
one_group = openmc.EnergyFilter([0., 1.e7])
dg_filter = openmc.DelayedGroupFilter(range(1,7))

# Create pin cell mesh
mesh = openmc.Mesh()
mesh.type = 'regular'
mesh.dimension = [51, 51]
mesh.lower_left  = [-32.13, -32.13]
mesh.upper_right = [ 32.13,  32.13]
pin_cell_mesh_filter = openmc.MeshFilter(mesh)

# Create a material filter
fuel_ids = []
mat_ids = []
for name,mat in materials.items():
    mat_ids.append(mat.id)
    if 'UO2' in name or 'MOX' in name:
        fuel_ids.append(mat.id)

fuel_filter = openmc.MaterialFilter(fuel_ids)
mat_filter = openmc.MaterialFilter(mat_ids)

# Pin powers tally
tallies['kappa-fission'] = openmc.Tally(name='kappa-fission')
tallies['kappa-fission'].scores = ['kappa-fission']
tallies['kappa-fission'].filters = [pin_cell_mesh_filter]

# Parameters for each material
tallies['inverse-velocity'] = openmc.Tally(name='inverse-velocity')
tallies['inverse-velocity'].scores = ['inverse-velocity', 'flux']
tallies['inverse-velocity'].filters = [mat_filter, seven_groups]

tallies['delayed-nu-fission'] = openmc.Tally(name='delayed-nu-fission')
tallies['delayed-nu-fission'].scores = ['delayed-nu-fission']
tallies['delayed-nu-fission'].filters = [fuel_filter, dg_filter]

tallies['nu-fission'] = openmc.Tally(name='nu-fission')
tallies['nu-fission'].scores = ['nu-fission']
tallies['nu-fission'].filters = [fuel_filter]

tallies['decay-rate'] = openmc.Tally(name='decay-rate')
tallies['decay-rate'].scores = ['decay-rate']
tallies['decay-rate'].filters = [fuel_filter, dg_filter]
