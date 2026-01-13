import matplotlib.pyplot as plt
import numpy as np
import openmc

###############################################################################
# Create materials for the problem

uo2 = openmc.Material(name="UO2 fuel at 2.4% wt enrichment")
uo2.set_density("g/cm3", 10.29769)
uo2.add_element("U", 1.0, enrichment=2.4)
uo2.add_element("O", 2.0)

helium = openmc.Material(name="Helium for gap")
helium.set_density("g/cm3", 0.001598)
helium.add_element("He", 2.4044e-4)

zircaloy = openmc.Material(name="Zircaloy 4")
zircaloy.set_density("g/cm3", 6.55)
zircaloy.add_element("Sn", 0.014, "wo")
zircaloy.add_element("Fe", 0.00165, "wo")
zircaloy.add_element("Cr", 0.001, "wo")
zircaloy.add_element("Zr", 0.98335, "wo")

borated_water = openmc.Material(name="Borated water")
borated_water.set_density("g/cm3", 0.740582)
borated_water.add_element("B", 2.0e-4)  # 3x the original pincell
borated_water.add_element("H", 5.0e-2)
borated_water.add_element("O", 2.4e-2)
borated_water.add_s_alpha_beta("c_H_in_H2O")

###############################################################################
# Define problem geometry

# Create cylindrical surfaces
fuel_or = openmc.ZCylinder(r=0.39218, name="Fuel OR")
clad_ir = openmc.ZCylinder(r=0.40005, name="Clad IR")
clad_or = openmc.ZCylinder(r=0.45720, name="Clad OR")

# Create a region represented as the inside of a rectangular prism
pitch = 1.25984
box = openmc.model.RectangularPrism(pitch, pitch, boundary_type="reflective")

# Create cells, mapping materials to regions
fuel = openmc.Cell(fill=uo2, region=-fuel_or)
gap = openmc.Cell(fill=helium, region=+fuel_or & -clad_ir)
clad = openmc.Cell(fill=zircaloy, region=+clad_ir & -clad_or)
water = openmc.Cell(fill=borated_water, region=+clad_or & -box)

# Create a model and assign geometry
model = openmc.Model()
model.geometry = openmc.Geometry([fuel, gap, clad, water])

###############################################################################
# Define problem settings

# Set the mode
model.settings.run_mode = "fixed source"

# Indicate how many batches and particles to run
model.settings.batches = 10
model.settings.particles = 10000

# Set time cutoff (we only care about t < 100 seconds, see tally below)
model.settings.cutoff = {"time_neutron": 100}

# Create the neutron pulse source (by default, isotropic direction, t=0)
space = openmc.stats.Point()  # At the origin (0, 0, 0)
energy = openmc.stats.delta_function(14.1e6)  # At 14.1 MeV
model.settings.source = openmc.IndependentSource(space=space, energy=energy)

###############################################################################
# Define tallies

# Create time filter
t_grid = np.insert(np.logspace(-6, 2, 100), 0, 0.0)
time_filter = openmc.TimeFilter(t_grid)

# Tally for total neutron density in time
density_tally = openmc.Tally(name="Density")
density_tally.filters = [time_filter]
density_tally.scores = ["inverse-velocity"]

# Add tallies to model
model.tallies = openmc.Tallies([density_tally])


# Run the model
model.run(apply_tally_results=True)

# Bin-averaged result
density_mean = density_tally.mean.ravel() / np.diff(t_grid)

# Plot particle density versus time
fig, ax = plt.subplots()
ax.stairs(density_mean, t_grid)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Total density")
ax.grid()
plt.show()
