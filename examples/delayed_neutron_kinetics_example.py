"""
Example demonstrating delayed neutron kinetics calculations in OpenMC.

This example shows how to enable calculation of:
- k_prompt: prompt neutron k-effective
- beta_eff: effective delayed neutron fraction
- alpha eigenvalues: reactor kinetics parameters (both k-based and rate-based)

The kinetics parameters are automatically calculated during an eigenvalue
simulation when enabled in the settings.
"""

import openmc

# Create a simple pin cell model for demonstration
# (In practice, use your own geometry/materials)

# Materials
uo2 = openmc.Material(name='UO2 fuel')
uo2.set_density('g/cm3', 10.0)
uo2.add_nuclide('U235', 0.03)
uo2.add_nuclide('U238', 0.97)
uo2.add_nuclide('O16', 2.0)

water = openmc.Material(name='Water')
water.set_density('g/cm3', 1.0)
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([uo2, water])
materials.export_to_xml()

# Geometry
fuel_radius = openmc.ZCylinder(r=0.39)
cell1 = openmc.Cell(fill=uo2, region=-fuel_radius)
cell2 = openmc.Cell(fill=water, region=+fuel_radius)
universe = openmc.Universe(cells=[cell1, cell2])
box = openmc.model.RectangularPrism(1.26, 1.26, boundary_type='reflective')
root_cell = openmc.Cell(fill=universe, region=-box)
geometry = openmc.Geometry([root_cell])
geometry.export_to_xml()

# Settings
settings = openmc.Settings()
settings.batches = 120
settings.inactive = 20
settings.particles = 10000

# Enable delayed neutron kinetics calculations
# This will calculate k_prompt and beta_eff
settings.calculate_prompt_k = True

# Enable alpha eigenvalue calculations
# This will also calculate alpha (k-based) and alpha (rate-based)
settings.calculate_alpha = True

settings.export_to_xml()

# Run the simulation
# openmc.run()

# Analyze results from statepoint
# sp = openmc.StatePoint('statepoint.120.h5')
# print(f"k-effective: {sp.keff}")
# print(f"k-prompt: {sp.k_prompt}")
# print(f"Beta-effective: {sp.beta_eff}")
# print(f"Prompt generation time: {sp.prompt_gen_time}")
# print(f"Alpha (k-based): {sp.alpha_k_based}")
# print(f"Alpha (rate-based): {sp.alpha_rate_based}")

print("Example files created successfully!")
print("To run this example:")
print("  1. Ensure you have cross section data available")
print("  2. Uncomment the openmc.run() line")
print("  3. Uncomment the analysis section")
print("  4. Run: python delayed_neutron_kinetics_example.py")
