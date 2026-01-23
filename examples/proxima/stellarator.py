import numpy as np
import openmc  # install with: conda install -c conda-forge openmc
import math

materials = openmc.Materials()

# simplified material definitions have been used to keen this example minimal
material_layer_1 = openmc.Material(name='layer_1')
material_layer_1.add_nuclide("Fe56", 1, "ao")
material_layer_1.set_density("g/cm3", 7.)
materials.append(material_layer_1)

material_layer_2 = openmc.Material(name='layer_2')
material_layer_2.add_nuclide("Li6", 0.9, "ao")
material_layer_2.add_nuclide("Li7", 0.1, "ao")
material_layer_2.set_density("g/cm3", 2.)
materials.append(material_layer_2)

material_layer_3 = openmc.Material(name='layer_3')
material_layer_3.add_nuclide("Fe56", 1, "ao")
material_layer_3.set_density("g/cm3", 7.)
materials.append(material_layer_3)

material_magnet = openmc.Material(name='magnet')
material_magnet.add_nuclide("Fe56", 1, "ao")
material_magnet.set_density("g/cm3", 7.)
materials.append(material_magnet)


bound_dag_univ = openmc.DAGMCUniverse(filename="dagmc_surface_mesh.h5m").bounded_universe()
geometry = openmc.Geometry(root=bound_dag_univ)

rmesh = openmc.RegularMesh.from_domain(bound_dag_univ)
rmesh.id = 2
rmesh.dimension = [200, 200, 200]
rmesh_filter = openmc.MeshFilter(rmesh)
rmesh_tally_h3 = openmc.Tally(name="H3_rmesh_tally")
rmesh_tally_h3.filters = [rmesh_filter]
rmesh_tally_h3.scores = ["H3-production"]

material_filter_layer_2 = openmc.MaterialFilter(material_layer_2)
material_tally_h3 = openmc.Tally(name="H3_material_tally")
material_tally_h3.filters = [material_filter_layer_2]
material_tally_h3.scores = ["H3-production"]

material_filter_magnet = openmc.MaterialFilter(material_magnet)
material_tally_magnet_heat = openmc.Tally(name="magnet_heating_material_tally")
material_tally_magnet_heat.filters = [material_filter_magnet]
material_tally_magnet_heat.scores = ["heating"]

# tallies = openmc.Tallies([material_tally_h3, material_tally_magnet_heat, mesh_tally_h3])
tallies = openmc.Tallies([material_tally_h3, material_tally_magnet_heat, rmesh_tally_h3])
# tallies = openmc.Tallies([mesh_tally, material_tally])


# initializes a new source object
my_source = openmc.IndependentSource()

# the distribution of radius is just a single value at the major radius
radius = openmc.stats.Discrete([1397.626385546529], [1])

# the distribution of source z values is just a single value
z_values = openmc.stats.Discrete([0], [1])

# the distribution of source azimuthal angles values is a uniform distribution between 0 and 2 Pi
angle = openmc.stats.Uniform(a=0., b=2* 3.14159265359)

# this makes the ring source using the three distributions and a radius
#my_source.space = openmc.stats.CylindricalIndependent(r=radius, phi=angle, z=z_values, origin=(0.0, 0.0, 0.0))

# sets the direction to isotropic
my_source.angle = openmc.stats.Isotropic()
# sets the energy distribution to 100% 14MeV neutrons, this is a simplification
# normally this would be a distribution of energies that can be made using packages like
# fusion_neutron_utils https://github.com/fusion-energy/fusion_neutron_utils
energy = openmc.stats.Discrete([14.06e6], [1])

# Let's make a ring of sources
sources = []
n_sources = 1000
#radius = 1397.626385546529
#radius = 1500.0
radius = 1480.0
z = 0.0

for i in range(n_sources):
    theta = 2.0 * math.pi * i / n_sources     # angle in radians
    x = radius * math.cos(theta)
    y = radius * math.sin(theta)

    src = openmc.IndependentSource()
    src.space = openmc.stats.Point([x, y, z])
    src.energy = energy
    sources.append(src)

# specifies the simulation computational intensity
settings = openmc.Settings()
settings.batches = 5  # this is set to 2 for testing purposes, this should be increased for an accurate simulation
settings.particles = 100  # this is set to 10 for testing purposes, this should be increased for an accurate simulation
settings.run_mode = "fixed source"
settings.source = sources
settings.photon_transport = True


# builds the openmc model
model = openmc.Model(
    materials=materials, geometry=geometry, settings=settings, tallies=tallies
)

#model.export_to_model_xml()

mesh = openmc.RegularMesh()
mesh.lower_left = geometry.bounding_box.lower_left
mesh.upper_right = geometry.bounding_box.upper_right
print("Lower Left:  ({}, {}, {})".format(mesh.lower_left[0], mesh.lower_left[1], mesh.lower_left[2]))
print("Upper Right: ({}, {}, {})".format(mesh.upper_right[0], mesh.upper_right[1], mesh.upper_right[2]))
resolution = 20.0
# Use the ceiling rather than the floor
n_x = math.ceil((mesh.upper_right[0] - mesh.lower_left[0]) / resolution)
n_y = math.ceil((mesh.upper_right[1] - mesh.lower_left[1]) / resolution)
n_z = math.ceil((mesh.upper_right[2] - mesh.lower_left[2]) / resolution)
mesh.dimension = [n_x, n_y, n_z]
print("Dimensions: ({}, {}, {})".format(n_x, n_y, n_z))

tally = openmc.Tally(name='ww_mesh_flux')
tally.filters = [openmc.MeshFilter(mesh)]
tally.scores = ['flux']
#model.tallies.append(tally)
model.tallies = [tally]

for mat in model.materials:
    mat.temperature = 300.0

model.convert_to_multigroup(
    method='stochastic_slab',
    groups='CASMO-4',
    nparticles=100000,
    overwrite_mgxs_library=False,
    correction=None
    #source_energy=14.6e6
)


model.settings.batches = 200
model.settings.inactive = 100
model.settings.particles = 30000
model.settings.random_ray['distance_inactive'] = 500.0
model.settings.random_ray['distance_active'] = 5000.0
model.settings.random_ray['volume_estimator'] = 'naive'
model.settings.random_ray['source_shape'] = 'flat'
#model.settings.random_ray['source_shape'] = 'linear'
model.settings.random_ray['volume_normalized_flux_tallies'] = False
settings.photon_transport = False

wwg = openmc.WeightWindowGenerator(
    mesh,
    method='fw_cadis',
    max_realizations=model.settings.batches
)

# Disable weight window generation for now
#model.settings.weight_window_generators = [wwg]


bbox      = geometry.bounding_box
print(bbox)
box_src = openmc.stats.Box(bbox.lower_left, bbox.upper_right)
rr_src  = openmc.IndependentSource(
    space       = box_src,
    constraints = {'domains': [geometry.root_universe]}
    )
model.settings.random_ray['ray_source'] = rr_src
model.settings.random_ray['source_region_meshes'] = [(mesh, [model.geometry.root_universe])]

model.geometry.remove_redundant_surfaces()

# Random Ray Plot
plot = openmc.Plot()
plot.origin = [0, 0, 0]
plot.width = [(geometry.bounding_box.upper_right[i] - geometry.bounding_box.lower_left[i]) for i in range(3)]
plot.pixels = [300, 300, 300]
plot.type = 'voxel'
model.plots = [plot]

model.export_to_model_xml()