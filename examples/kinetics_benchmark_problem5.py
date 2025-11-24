"""
Reactor Kinetics Benchmark Problem 5: Sub-critical Fast Neutron System

This benchmark problem calculates delayed neutron kinetics parameters for a
sub-critical uranium sphere with exactly half the density of Problem 1 (Godiva).
This problem demonstrates the difference between static and dynamic criticality
calculations when the system is significantly sub-critical.

The problem demonstrates:
- Calculation of k_eff and k_prompt for a deeply sub-critical system
- How beta_eff behaves in sub-critical configurations
- Significant differences between static and dynamic alpha eigenvalues
- Importance of prompt vs delayed neutron contributions in sub-critical systems
- When the system is far from critical, dynamic results are physically appropriate

Reference:
Cullen, Dermott E., Christopher J. Clouse, Richard Procassini, and Robert C. Little.
Static and Dynamic Criticality: Are They Different? UCRL-TR-201506.
Livermore, CA: Lawrence Livermore National Laboratory, November 22, 2003.

Author: William Zywiec (willzywiec@gmail.com)
"""

import openmc
import numpy as np

print("=" * 70)
print("Problem 5: Sub-critical Fast Neutron System")
print("=" * 70)

# =============================================================================
# Materials
# =============================================================================

# HEU sphere at half the density of Problem 1 - exact specification from Cullen et al. (2003)
# Density: 9.3999 grams/cc (exactly half of Problem 1: 18.7398 / 2)
# Composition (atom fractions - same as Problem 1):
#   U-235: 0.937695
#   U-238: 0.052053
#   U-234: 0.010252
heu_subcritical = openmc.Material(name='HEU Sub-critical')
heu_subcritical.set_density('g/cm3', 9.3999)
heu_subcritical.add_nuclide('U235', 0.937695)
heu_subcritical.add_nuclide('U238', 0.052053)
heu_subcritical.add_nuclide('U234', 0.010252)
heu_subcritical.temperature = 293.6  # Room temperature

materials = openmc.Materials([heu_subcritical])
materials.export_to_xml()

print("\nMaterial: HEU at Half Density (UCRL-TR-201506)")
print(f"  Density: 9.3999 g/cm³ (half of Problem 1)")
print(f"  U-235: 0.937695 (atom fraction)")
print(f"  U-238: 0.052053 (atom fraction)")
print(f"  U-234: 0.010252 (atom fraction)")
print(f"  Temperature: 293.6 K")
print(f"  Note: Same isotopic composition as Problem 1")

# =============================================================================
# Geometry
# =============================================================================

# Same geometry as Problem 1 - only density differs
radius = 8.7407  # cm (same as Problem 1)

sphere_surface = openmc.Sphere(r=radius, boundary_type='vacuum')
sphere_cell = openmc.Cell(fill=heu_subcritical, region=-sphere_surface)

geometry = openmc.Geometry([sphere_cell])
geometry.export_to_xml()

print(f"\nGeometry: Bare sphere with vacuum boundary")
print(f"  Radius: {radius} cm (same as Problem 1)")
print(f"  Characteristics: Homogeneous, fast neutron, very sub-critical")
print(f"  This problem demonstrates static vs dynamic criticality")
print(f"  when the system is far from critical.")

# =============================================================================
# Settings
# =============================================================================

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.particles = 50000
settings.batches = 150
settings.inactive = 50

# Use a point source at the center of the sphere
settings.source = openmc.IndependentSource(space=openmc.stats.Point((0, 0, 0)))

# Enable delayed neutron kinetics calculations
settings.calculate_prompt_k = True
settings.calculate_alpha = True

settings.export_to_xml()

print(f"\nSimulation Settings:")
print(f"  Particles per batch: {settings.particles}")
print(f"  Total batches: {settings.batches}")
print(f"  Inactive batches: {settings.inactive}")
print(f"  Kinetics calculations: ENABLED")
print(f"")
print(f"Note: For sub-critical systems, convergence may require more particles")

# =============================================================================
# Run Simulation
# =============================================================================

print("\n" + "=" * 70)
print("Running simulation...")
print("=" * 70)

# Uncomment to run the simulation
openmc.run()

