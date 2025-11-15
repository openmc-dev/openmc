"""
Reactor Kinetics Benchmark Problem 1: Godiva (Near Critical, Fast Neutron System)

This benchmark problem calculates delayed neutron kinetics parameters for the
Godiva bare highly enriched uranium sphere - a classic criticality benchmark
that is near prompt critical. This simple geometry allows for validation of
k_prompt, beta_eff, and alpha eigenvalue calculations.

The problem demonstrates:
- Calculation of k_eff and k_prompt for a near-critical fast system
- Determination of beta_eff from the difference between k_eff and k_prompt
- Alpha eigenvalue calculations using both k-based and rate-based methods
- Comparison of static vs dynamic criticality approaches
- Validation that beta_eff matches known delayed neutron fraction data

Reference:
Cullen, Dermott E., Christopher J. Clouse, Richard Procassini, and Robert C. Little.
Static and Dynamic Criticality: Are They Different? UCRL-TR-201506.
Livermore, CA: Lawrence Livermore National Laboratory, November 22, 2003.

Author: William Zywiec (willzywiec@gmail.com)
"""

import openmc
import numpy as np

print("=" * 70)
print("Problem 1: Godiva (Near Critical, Fast Neutron System)")
print("=" * 70)

# =============================================================================
# Materials
# =============================================================================

# Godiva HEU sphere - exact specification from Cullen et al. (2003)
# Density: 18.7398 grams/cc
# Composition (atom fractions):
#   U-235: 0.937695
#   U-238: 0.052053
#   U-234: 0.010252
godiva = openmc.Material(name='Godiva HEU')
godiva.set_density('g/cm3', 18.7398)
godiva.add_nuclide('U235', 0.937695)
godiva.add_nuclide('U238', 0.052053)
godiva.add_nuclide('U234', 0.010252)
godiva.temperature = 293.6  # Room temperature

materials = openmc.Materials([godiva])
materials.export_to_xml()

print("\nMaterial: Godiva HEU (UCRL-TR-201506)")
print(f"  Density: 18.7398 g/cm³")
print(f"  U-235: 0.937695 (atom fraction)")
print(f"  U-238: 0.052053 (atom fraction)")
print(f"  U-234: 0.010252 (atom fraction)")
print(f"  Temperature: 293.6 K")

# =============================================================================
# Geometry
# =============================================================================

# Godiva sphere radius from Cullen et al. (2003)
radius = 8.7407  # cm

sphere_surface = openmc.Sphere(r=radius, boundary_type='vacuum')
sphere_cell = openmc.Cell(fill=godiva, region=-sphere_surface)

geometry = openmc.Geometry([sphere_cell])
geometry.export_to_xml()

print(f"\nGeometry: Bare sphere with vacuum boundary")
print(f"  Radius: {radius} cm (UCRL-TR-201506)")
print(f"  Characteristics: Homogeneous, fast neutron, near prompt critical")

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

# =============================================================================
# Run Simulation
# =============================================================================

print("\n" + "=" * 70)
print("Running simulation...")
print("=" * 70)

# Uncomment to run the simulation
# openmc.run()

# =============================================================================
# Analysis
# =============================================================================

print("\nTo analyze results after running, use:")
print("-" * 70)
print("""
import openmc

sp = openmc.StatePoint('statepoint.150.h5')

print("\\nResults:")
print("=" * 70)
print(f"k-effective:              {sp.keff}")
print(f"k-prompt:                 {sp.k_prompt}")
print(f"Beta-effective:           {sp.beta_eff}")
print(f"Prompt generation time:   {sp.prompt_gen_time} seconds")
print(f"Alpha (k-based):          {sp.alpha_k_based} 1/s")
print(f"Alpha (rate-based):       {sp.alpha_rate_based} 1/s")

# Expected results for Godiva (from Cullen et al. 2003):
# - k_eff should be close to 1.0 (near critical)
# - k_prompt should be less than k_eff
# - beta_eff should be approximately 0.0065-0.0070 (0.65-0.70%)
# - alpha should be small (near zero) for a near-critical system
# - Both alpha methods should agree within uncertainties

print("\\nExpected values for Godiva (UCRL-TR-201506):")
print("  k_eff ≈ 1.0 (near critical)")
print("  beta_eff ≈ 0.0065-0.0070 (0.65-0.70%)")
print("  k_prompt ≈ 0.993-0.994")
print("  This problem demonstrates static vs dynamic criticality concepts.")
""")
print("=" * 70)

print("\nBenchmark Problem 1 files created successfully!")
print("To run: python kinetics_benchmark_problem1.py")
print("(Requires ENDF/B cross section data)")
