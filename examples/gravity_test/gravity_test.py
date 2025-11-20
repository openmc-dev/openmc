"""
Simple example demonstrating gravitational acceleration in OpenMC.

This example creates a simple geometry with a neutron source and enables
gravity to observe its effect on particle trajectories.
"""

import openmc

# Create materials
mat = openmc.Material()
mat.add_nuclide('H1', 1.0)
mat.set_density('g/cm3', 0.001)  # Very low density for long free paths
materials = openmc.Materials([mat])
materials.export_to_xml()

# Create geometry - simple box
box = openmc.Cell()
box.region = -openmc.XPlane(100) & +openmc.XPlane(-100) & \
             -openmc.YPlane(100) & +openmc.YPlane(-100) & \
             -openmc.ZPlane(100) & +openmc.ZPlane(-100)
box.fill = mat

root = openmc.Universe(cells=[box])
geometry = openmc.Geometry(root)
geometry.export_to_xml()

# Create source at the top, shooting horizontally
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0, 0, 50))  # Start at z=50
source.angle = openmc.stats.Monodirectional((1, 0, 0))  # Shoot in +x direction
source.energy = openmc.stats.Discrete([1.0e6], [1.0])  # 1 MeV neutron

# Settings
settings = openmc.Settings()
settings.batches = 10
settings.inactive = 0
settings.particles = 100
settings.source = source
settings.run_mode = 'fixed source'

# Enable gravity with strong field for demonstration
# Use 1000x Earth gravity to make effects visible
settings.gravity = {
    'enabled': True,
    'acceleration': [0.0, 0.0, -980000.0]  # Strong downward gravity (cm/s^2)
}

settings.export_to_xml()

# Create tallies to track particle positions
mesh = openmc.RegularMesh()
mesh.lower_left = (-100, -100, -100)
mesh.upper_right = (100, 100, 100)
mesh.dimension = (20, 1, 20)  # Fine mesh in x-z plane

mesh_filter = openmc.MeshFilter(mesh)
tally = openmc.Tally()
tally.filters = [mesh_filter]
tally.scores = ['flux']

tallies = openmc.Tallies([tally])
tallies.export_to_xml()

print("Gravity test setup complete!")
print("Run with: openmc")
print("")
print("Expected behavior:")
print("- Neutrons start at z=50, traveling in +x direction")
print("- With gravity in -z direction, particles should curve downward")
print("- The flux mesh will show the curved trajectory")
print("")
print("To compare, run once with gravity enabled and once with:")
print("  settings.gravity['enabled'] = False")
