from math import pi
import os
import openmc
import openmc.deplete
import openmc.data
import matplotlib.pyplot as plt

openmc.config['cross_sections'] = "/mnt/d/OpenMC_debug/evaluations/endfb-viii.0-hdf5/cross_sections.xml"

#initiate an instance of the Chain class
mychain = openmc.deplete.Chain()

#Find all relevant endf files
decay_files = []
decay_folder = "/home/shauksson/evaluations/ENDF-B-VIII.0/decay/"
for file in os.listdir(decay_folder):
    if file[-4:] == "endf":
        decay_files.append(decay_folder + file)


nfy_files = []
nfy_folder = "/home/shauksson/evaluations/ENDF-B-VIII.0/nfy/"
for file in os.listdir(nfy_folder):
    if file[-4:] == "endf":
        nfy_files.append(nfy_folder + file)


xs_files = []
xs_folder = "/home/shauksson/evaluations/ENDF-B-VIII.0/neutrons/"
for file in os.listdir(xs_folder):
    if file[-4:] == "endf":
        xs_files.append(xs_folder + file)

sfy_files = []
sfy_folder = "/home/shauksson/evaluations/ENDF-B-VIII.0/sfy/"
for file in os.listdir(sfy_folder):
    if file[-4:] == "endf":
        sfy_files.append(sfy_folder + file)

#read from endf file
mychain = mychain.from_endf(decay_files,nfy_files,xs_files,sfy_files)

#Print chain to xml file
chain_file = "./chain.xml"
mychain.export_to_xml(chain_file)


###############################################################################
#                              Define materials
###############################################################################

# Instantiate some Materials and register the appropriate Nuclides
uo2 = openmc.Material(name='UO2 fuel at 2.4% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)

helium = openmc.Material(name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

zircaloy = openmc.Material(name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014, 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001, 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')

borated_water = openmc.Material(name='Borated water')
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_element('B', 4.0e-5)
borated_water.add_element('H', 5.0e-2)
borated_water.add_element('O', 2.4e-2)
borated_water.add_s_alpha_beta('c_H_in_H2O')

###############################################################################
#                             Create geometry
###############################################################################

# Define surfaces
pitch = 1.25984
fuel_or = openmc.ZCylinder(r=0.39218, name='Fuel OR')
clad_ir = openmc.ZCylinder(r=0.40005, name='Clad IR')
clad_or = openmc.ZCylinder(r=0.45720, name='Clad OR')
box = openmc.model.RectangularPrism(pitch, pitch, boundary_type='reflective')

# Define cells
fuel = openmc.Cell(fill=uo2, region=-fuel_or)
gap = openmc.Cell(fill=helium, region=+fuel_or & -clad_ir)
clad = openmc.Cell(fill=zircaloy, region=+clad_ir & -clad_or)
water = openmc.Cell(fill=borated_water, region=+clad_or & -box)

# Define overall geometry
geometry = openmc.Geometry([fuel, gap, clad, water])

###############################################################################
#                     Set volumes of depletable materials
###############################################################################

# Set material volume for depletion. For 2D simulations, this should be an area.
uo2.volume = pi * fuel_or.r**2

###############################################################################
#                     Transport calculation settings
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings = openmc.Settings()
settings.batches = 20
settings.inactive = 10
settings.particles = 1000

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.IndependentSource(
    space=uniform_dist, constraints={'fissionable': True})

entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = [-0.39218, -0.39218, -1.e50]
entropy_mesh.upper_right = [0.39218, 0.39218, 1.e50]
entropy_mesh.dimension = [10, 10, 1]
settings.entropy_mesh = entropy_mesh

###########################

#Create model
model = openmc.Model(geometry=geometry, settings=settings)

# Create depletion "operator"
op = openmc.deplete.CoupledOperator(model, chain_file)

# Perform simulation using the predictor algorithm
time_steps = [1.0, 1.0, 1.0, 1.0, 1.0]  # days
power = 174  # W/cm, for 2D simulations only (use W for 3D)
integrator = openmc.deplete.PredictorIntegrator(op, time_steps, power, timestep_units='d')
integrator.integrate()
