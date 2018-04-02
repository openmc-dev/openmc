import openmc
import openmc.deplete
import numpy as np

###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 100
inactive = 10
particles = 1000

# Depletion simulation parameters
time_step = 1*24*60*60 # s
final_time = 15*24*60*60 # s
time_steps = np.full(final_time // time_step, time_step)

chain_file = './chain_simple.xml'
power = 174 # W/cm, for 2D simulations only (use W for 3D)

###############################################################################
#                              Define materials
###############################################################################

# Instantiate some Materials and register the appropriate Nuclides
uo2 = openmc.Material(material_id=1, name='UO2 fuel at 2.4% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)
uo2.depletable = True

helium = openmc.Material(material_id=2, name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

zircaloy = openmc.Material(material_id=3, name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014  , 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001  , 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')

borated_water = openmc.Material(material_id=4, name='Borated water')
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_element('B', 4.0e-5)
borated_water.add_element('H', 5.0e-2)
borated_water.add_element('O', 2.4e-2)
borated_water.add_s_alpha_beta('c_H_in_H2O')

###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate ZCylinder surfaces
fuel_or = openmc.ZCylinder(surface_id=1, x0=0, y0=0, R=0.39218, name='Fuel OR')
clad_ir = openmc.ZCylinder(surface_id=2, x0=0, y0=0, R=0.40005, name='Clad IR')
clad_or = openmc.ZCylinder(surface_id=3, x0=0, y0=0, R=0.45720, name='Clad OR')
left = openmc.XPlane(surface_id=4, x0=-0.62992, name='left')
right = openmc.XPlane(surface_id=5, x0=0.62992, name='right')
bottom = openmc.YPlane(surface_id=6, y0=-0.62992, name='bottom')
top = openmc.YPlane(surface_id=7, y0=0.62992, name='top')

left.boundary_type = 'reflective'
right.boundary_type = 'reflective'
top.boundary_type = 'reflective'
bottom.boundary_type = 'reflective'

# Instantiate Cells
fuel = openmc.Cell(cell_id=1, name='cell 1')
gap = openmc.Cell(cell_id=2, name='cell 2')
clad = openmc.Cell(cell_id=3, name='cell 3')
water = openmc.Cell(cell_id=4, name='cell 4')

# Use surface half-spaces to define regions
fuel.region = -fuel_or
gap.region = +fuel_or & -clad_ir
clad.region = +clad_ir & -clad_or
water.region = +clad_or & +left & -right & +bottom & -top

# Register Materials with Cells
fuel.fill = uo2
gap.fill = helium
clad.fill = zircaloy
water.fill = borated_water

# Instantiate Universe
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
root.add_cells([fuel, gap, clad, water])

# Instantiate a Geometry, register the root Universe
geometry = openmc.Geometry(root)

###############################################################################
#                   Exporting to OpenMC materials.xml file
###############################################################################

# Compute cell areas
area = {}
area[fuel] = np.pi * fuel_or.coefficients['R'] ** 2

# Set materials volume for depletion. Set to an area for 2D simulations
uo2.volume = area[fuel]

###############################################################################
#                   Exporting to OpenMC settings.xml file
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
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

op = openmc.deplete.Operator(geometry, settings_file, chain_file)

# Perform simulation using the predictor algorithm
openmc.deplete.integrator.predictor(op, time_steps, power)

###############################################################################
#                    Read depletion calculation results
###############################################################################

# Open results file
results = openmc.deplete.ResultsList("depletion_results.h5")

# Obtain K_eff as a function of time
time, keff = results.get_eigenvalue()
                                                                                                                                                                                                                                                                                                                                                                                                                                                        
# Obtain U235 concentration as a function of time
time, n_U235 = results.get_atoms('1', 'U235')
