from math import pi

import openmc
import openmc.deplete
import matplotlib.pyplot as plt

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
box = openmc.model.rectangular_prism(pitch, pitch, boundary_type='reflective')

# Define cells
fuel = openmc.Cell(fill=uo2, region=-fuel_or)
gap = openmc.Cell(fill=helium, region=+fuel_or & -clad_ir)
clad = openmc.Cell(fill=zircaloy, region=+clad_ir & -clad_or)
water = openmc.Cell(fill=borated_water, region=+clad_or & box)

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
settings.batches = 100
settings.inactive = 10
settings.particles = 1000

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)

entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = [-0.39218, -0.39218, -1.e50]
entropy_mesh.upper_right = [0.39218, 0.39218, 1.e50]
entropy_mesh.dimension = [10, 10, 1]
settings.entropy_mesh = entropy_mesh

###############################################################################
#                   Initialize and run depletion calculation
###############################################################################

# Create depletion "operator"
chain_file = './chain_simple.xml'
op = openmc.deplete.Operator(geometry, settings, chain_file)

# Perform simulation using the predictor algorithm
time_steps = [1.0, 1.0, 1.0, 1.0, 1.0]  # days
power = 174  # W/cm, for 2D simulations only (use W for 3D)
integrator = openmc.deplete.PredictorIntegrator(op, time_steps, power, timestep_units='d')
integrator.integrate()

###############################################################################
#                    Read depletion calculation results
###############################################################################

# Open results file
results = openmc.deplete.ResultsList.from_hdf5("depletion_results.h5")

# Obtain K_eff as a function of time
time, keff = results.get_eigenvalue()

# Obtain U235 concentration as a function of time
time, n_U235 = results.get_atoms('1', 'U235')

# Obtain Xe135 capture reaction rate as a function of time
time, Xe_capture = results.get_reaction_rate('1', 'Xe135', '(n,gamma)')

###############################################################################
#                            Generate plots
###############################################################################

days = 24*60*60
plt.figure()
plt.plot(time/days, keff, label="K-effective")
plt.xlabel("Time (days)")
plt.ylabel("Keff")
plt.show()

plt.figure()
plt.plot(time/days, n_U235, label="U235")
plt.xlabel("Time (days)")
plt.ylabel("n U5 (-)")
plt.show()

plt.figure()
plt.plot(time/days, Xe_capture, label="Xe135 capture")
plt.xlabel("Time (days)")
plt.ylabel("RR (-)")
plt.show()
plt.close('all')
