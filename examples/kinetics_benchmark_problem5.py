"""
Reactor Kinetics Benchmark Problem 5: PWR Fuel Pin Cell

This benchmark problem calculates delayed neutron kinetics parameters for a
typical Pressurized Water Reactor (PWR) fuel pin cell with water moderator.
This represents a more realistic reactor configuration with thermal neutron
spectrum and strong moderation effects.

The problem demonstrates:
- Kinetics parameters in a thermal spectrum environment
- Effect of moderation on delayed neutron behavior
- Comparison of beta_eff between fast and thermal systems
- Alpha eigenvalue sensitivity to spectrum

Reference: Based on typical PWR pin cell benchmarks
Author: William Zywiec (willzywiec@gmail.com)
"""

import openmc
import numpy as np

print("=" * 70)
print("Reactor Kinetics Benchmark Problem 5: PWR Fuel Pin Cell")
print("=" * 70)

# =============================================================================
# Materials
# =============================================================================

# UO2 fuel at 3.1% enrichment (typical PWR fresh fuel)
uo2 = openmc.Material(name='UO2 Fuel 3.1%')
uo2.set_density('g/cm3', 10.29)
uo2.add_nuclide('U235', 0.031)
uo2.add_nuclide('U238', 0.969)
uo2.add_element('O', 2.0)
uo2.temperature = 900.0  # Fuel temperature in K

# Zircaloy-4 cladding
zirconium = openmc.Material(name='Zircaloy-4')
zirconium.set_density('g/cm3', 6.56)
zirconium.add_element('Zr', 0.9821)
zirconium.add_element('Sn', 0.0145)
zirconium.add_element('Fe', 0.0021)
zirconium.add_element('Cr', 0.0010)
zirconium.add_element('O', 0.0012)
zirconium.temperature = 600.0

# Water moderator with boron (typical PWR conditions)
water = openmc.Material(name='Light Water with Boron')
water.set_density('g/cm3', 0.7405)  # Density at PWR conditions
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.add_element('B', 650e-6)  # 650 ppm boron
water.add_s_alpha_beta('c_H_in_H2O')
water.temperature = 600.0

materials = openmc.Materials([uo2, zirconium, water])
materials.export_to_xml()

print("\nMaterials:")
print("  UO2 Fuel:")
print("    Enrichment: 3.1% U-235")
print("    Density: 10.29 g/cm³")
print("    Temperature: 900 K")
print("  Zircaloy-4 Cladding:")
print("    Density: 6.56 g/cm³")
print("    Temperature: 600 K")
print("  Water Moderator:")
print("    Density: 0.7405 g/cm³")
print("    Boron concentration: 650 ppm")
print("    Temperature: 600 K")

# =============================================================================
# Geometry
# =============================================================================

# Typical PWR 17x17 lattice dimensions
fuel_or = 0.4096   # Fuel outer radius (cm)
clad_ir = 0.4178   # Cladding inner radius (cm)
clad_or = 0.4750   # Cladding outer radius (cm)
pitch = 1.26       # Pin pitch (cm)

# Surfaces
fuel_surface = openmc.ZCylinder(r=fuel_or)
clad_inner_surface = openmc.ZCylinder(r=clad_ir)
clad_outer_surface = openmc.ZCylinder(r=clad_or)

# Cells
fuel_cell = openmc.Cell(fill=uo2, region=-fuel_surface, name='Fuel')
gap_cell = openmc.Cell(region=+fuel_surface & -clad_inner_surface, name='Gap')
clad_cell = openmc.Cell(fill=zirconium, region=+clad_inner_surface & -clad_outer_surface, name='Cladding')
moderator_cell = openmc.Cell(fill=water, region=+clad_outer_surface, name='Moderator')

# Pin cell universe
pin_universe = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, moderator_cell])

# Create infinite lattice with reflective boundaries for pin cell calculation
boundary_box = openmc.model.RectangularPrism(
    width=pitch,
    height=pitch,
    boundary_type='reflective'
)
root_cell = openmc.Cell(fill=pin_universe, region=-boundary_box)

geometry = openmc.Geometry([root_cell])
geometry.export_to_xml()

print(f"\nGeometry: PWR Pin Cell (infinite lattice)")
print(f"  Fuel pellet radius: {fuel_or} cm")
print(f"  Cladding inner radius: {clad_ir} cm")
print(f"  Cladding outer radius: {clad_or} cm")
print(f"  Pin pitch: {pitch} cm")
print(f"  Boundary: Reflective (infinite lattice)")

# =============================================================================
# Settings
# =============================================================================

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.particles = 10000
settings.batches = 150
settings.inactive = 30

# Distributed source throughout the fuel
lower_left = (-pitch/2, -pitch/2, -1)
upper_right = (pitch/2, pitch/2, 1)
uniform_dist = openmc.stats.Box(lower_left, upper_right)
settings.source = openmc.IndependentSource(space=uniform_dist)

# Enable delayed neutron kinetics calculations
settings.calculate_prompt_k = True
settings.calculate_alpha = True

# Use Shannon entropy to monitor source convergence
entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = (-pitch/2, -pitch/2, -1)
entropy_mesh.upper_right = (pitch/2, pitch/2, 1)
entropy_mesh.dimension = (5, 5, 1)
settings.entropy_mesh = entropy_mesh

settings.export_to_xml()

print(f"\nSimulation Settings:")
print(f"  Particles per batch: {settings.particles}")
print(f"  Total batches: {settings.batches}")
print(f"  Inactive batches: {settings.inactive}")
print(f"  Kinetics calculations: ENABLED")
print(f"  Shannon entropy: ENABLED")

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
import matplotlib.pyplot as plt

sp = openmc.StatePoint('statepoint.150.h5')

print("\\nResults:")
print("=" * 70)
print(f"k-effective:              {sp.keff}")
print(f"k-prompt:                 {sp.k_prompt}")
print(f"Beta-effective:           {sp.beta_eff}")
print(f"Prompt generation time:   {sp.prompt_gen_time} seconds")
print(f"Alpha (k-based):          {sp.alpha_k_based} 1/s")
print(f"Alpha (rate-based):       {sp.alpha_rate_based} 1/s")

# Expected results for PWR systems:
# - k_inf typically 1.2-1.4 for fresh fuel
# - beta_eff approximately 0.0065-0.0075 (similar to fast systems)
# - Prompt generation time much longer than fast systems (~50-100 μs)
# - Alpha more negative due to longer generation time

print("\\nExpected values for PWR pin cell:")
print("  k_inf ≈ 1.2-1.4 (fresh fuel)")
print("  beta_eff ≈ 0.0065-0.0075 (0.65-0.75%)")
print("  Generation time ≈ 50-100 microseconds")
print("  Alpha will depend on k_inf and generation time")

# Plot k_eff and k_prompt convergence
plt.figure(figsize=(10, 6))
generations = range(len(sp.k_generation))
plt.plot(generations, sp.k_generation, label='k-effective', alpha=0.7)
plt.plot(generations, sp.k_prompt_generation, label='k-prompt', alpha=0.7)
plt.axhline(y=float(sp.keff.nominal_value), color='b', linestyle='--',
            label=f'Mean k-eff = {sp.keff}')
plt.axhline(y=float(sp.k_prompt.nominal_value), color='r', linestyle='--',
            label=f'Mean k-prompt = {sp.k_prompt}')
plt.xlabel('Generation')
plt.ylabel('k-eigenvalue')
plt.title('PWR Pin Cell: k-effective and k-prompt Convergence')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('pwr_kinetics_convergence.png', dpi=300)
print("\\nPlot saved: pwr_kinetics_convergence.png")

# Calculate beta_eff directly from difference
beta_calculated = (sp.keff - sp.k_prompt) / sp.keff
print(f"\\nBeta-eff (from k difference): {beta_calculated}")
print(f"Beta-eff (from simulation):   {sp.beta_eff}")
print("These should match within statistical uncertainties.")
""")
print("=" * 70)

print("\nBenchmark Problem 5 files created successfully!")
print("To run: python kinetics_benchmark_problem5.py")
print("(Requires ENDF/B cross section data with thermal scattering)")
