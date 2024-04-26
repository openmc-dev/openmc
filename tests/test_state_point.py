import openmc

# Define a material.
fuel = openmc.Material(name='UO2 Fuel')
fuel.add_element('U', 1.0, enrichment=4.0)
fuel.add_element('O', 2.0)
fuel.set_density('g/cm3', 10.0)

# Add the material to a materials collection.
materials = openmc.Materials([fuel])
materials.export_to_xml()

# Define surfaces.
fuel_or = openmc.ZCylinder(r=0.4, name='Fuel OR')
boundary = openmc.model.rectangular_prism(width=2.0, height=2.0, boundary_type='reflective')

# Define cells.
fuel_cell = openmc.Cell(name='fuel cell', fill=fuel, region=-fuel_or)
outside_fuel_cell = openmc.Cell(name='outside fuel', region=+fuel_or & -boundary)

# Add cells to universe.
universe = openmc.Universe(cells=[fuel_cell, outside_fuel_cell])
geometry = openmc.Geometry(universe)

geometry.export_to_xml()

# Create a tally to measure neutron flux.
tally = openmc.Tally(name='flux')
tally.scores = ['flux']
tallies = openmc.Tallies([tally])
tallies.export_to_xml()

# Define simulation settings.
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = 1000
settings.source = openmc.Source(space=openmc.stats.Point((0, 0, 0)))

settings.export_to_xml()

# Run the simulation
openmc.run()