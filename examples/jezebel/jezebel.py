import openmc

# Create plutonium metal material
pu = openmc.Material()
pu.set_density('sum')
pu.add_nuclide('Pu239', 3.7047e-02)
pu.add_nuclide('Pu240', 1.7512e-03)
pu.add_nuclide('Pu241', 1.1674e-04)
pu.add_element('Ga', 1.3752e-03)
mats = openmc.Materials([pu])
mats.export_to_xml()

# Create a single cell filled with the Pu metal
sphere = openmc.Sphere(r=6.3849, boundary_type='vacuum')
cell = openmc.Cell(fill=pu, region=-sphere)
geom = openmc.Geometry([cell])
geom.export_to_xml()

# Finally, define some run settings
settings = openmc.Settings()
settings.batches = 200
settings.inactive = 10
settings.particles = 10000
settings.export_to_xml()

# Run the simulation
openmc.run()

# Get the resulting k-effective value
n = settings.batches
with openmc.StatePoint(f'statepoint.{n}.h5') as sp:
    keff = sp.k_combined
    print(f'Final k-effective = {keff}')
