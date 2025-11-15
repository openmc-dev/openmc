"""
Reactor Kinetics Benchmark Problem 1: Bare Sphere Critical Assembly

This benchmark problem calculates delayed neutron kinetics parameters for a
bare homogeneous sphere of fissile material. This simple geometry allows for
validation of k_prompt, beta_eff, and alpha eigenvalue calculations against
analytical or reference Monte Carlo solutions.

The problem demonstrates:
- Calculation of k_eff and k_prompt for a critical system
- Determination of beta_eff from the difference between k_eff and k_prompt
- Alpha eigenvalue calculations using both k-based and rate-based methods
- Validation that beta_eff matches known delayed neutron fraction data

Reference: Based on typical bare sphere criticality benchmarks
Author: William Zywiec (willzywiec@gmail.com)
"""

import openmc
import numpy as np

print("=" * 70)
print("Reactor Kinetics Benchmark Problem 1: Bare Sphere Critical Assembly")
print("=" * 70)

# =============================================================================
# Materials
# =============================================================================

# Highly enriched uranium metal sphere (similar to Godiva)
# Density: 18.75 g/cm³, 93.71% U-235 enrichment
heu = openmc.Material(name='HEU Metal')
heu.set_density('g/cm3', 18.75)
heu.add_nuclide('U235', 0.9371)
heu.add_nuclide('U238', 0.0629)
heu.temperature = 293.6  # Room temperature

materials = openmc.Materials([heu])
materials.export_to_xml()

print("\nMaterial: HEU Metal Sphere")
print(f"  Density: 18.75 g/cm³")
print(f"  U-235 enrichment: 93.71%")
print(f"  Temperature: 293.6 K")

# =============================================================================
# Geometry
# =============================================================================

# Critical radius for bare HEU sphere is approximately 8.7 cm
# Adjust radius to achieve criticality (k_eff ≈ 1.0)
critical_radius = 8.71  # cm

sphere_surface = openmc.Sphere(r=critical_radius, boundary_type='vacuum')
sphere_cell = openmc.Cell(fill=heu, region=-sphere_surface)

geometry = openmc.Geometry([sphere_cell])
geometry.export_to_xml()

print(f"\nGeometry: Bare sphere with vacuum boundary")
print(f"  Critical radius: {critical_radius} cm")

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

# Expected results for HEU systems:
# - k_eff should be close to 1.0 (critical)
# - k_prompt should be less than k_eff
# - beta_eff should be approximately 0.0065-0.0070 (0.65-0.70%)
# - alpha should be negative for a critical system
# - Both alpha methods should agree within uncertainties

print("\\nExpected values for HEU:")
print("  beta_eff ≈ 0.0065-0.0070 (0.65-0.70%)")
print("  k_prompt ≈ 0.993-0.994 for k_eff = 1.0")
""")
print("=" * 70)

print("\nBenchmark Problem 1 files created successfully!")
print("To run: python kinetics_benchmark_problem1.py")
print("(Requires ENDF/B cross section data)")
