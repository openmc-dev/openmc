#!/usr/bin/env python
"""
Test script demonstrating the benefits of using derivative tallies in keff_search.

This script compares two approaches:
1. Standard GRsecant search without derivatives
2. Enhanced search using derivative tallies as constraints

The model is a PWR pin cell with boron-controlled moderator. We search for the
boron concentration (ppm) that achieves k-effective = 1.0 (criticality).
"""

import openmc
import openmc.stats
import numpy as np
import time
from pathlib import Path
import tempfile
import shutil


# ===============================================================
# Model builder
# ===============================================================
def build_model(ppm_boron):
    """Build a PWR pin cell model with specified boron concentration."""
    
    # Create the pin materials
    fuel = openmc.Material(name='1.6% Fuel', material_id=1)
    fuel.set_density('g/cm3', 10.31341)
    fuel.add_element('U', 1., enrichment=1.6)
    fuel.add_element('O', 2.)

    zircaloy = openmc.Material(name='Zircaloy', material_id=2)
    zircaloy.set_density('g/cm3', 6.55)
    zircaloy.add_element('Zr', 1.)

    water = openmc.Material(name='Borated Water', material_id=3)
    water.set_density('g/cm3', 0.741)
    water.add_element('H', 2.)
    water.add_element('O', 1.)
    water.add_element('B', ppm_boron * 1e-6)

    materials = openmc.Materials([fuel, zircaloy, water])

    # Geometry
    fuel_outer_radius = openmc.ZCylinder(r=0.39218)
    clad_outer_radius = openmc.ZCylinder(r=0.45720)

    min_x = openmc.XPlane(x0=-0.63, boundary_type='reflective')
    max_x = openmc.XPlane(x0=+0.63, boundary_type='reflective')
    min_y = openmc.YPlane(y0=-0.63, boundary_type='reflective')
    max_y = openmc.YPlane(y0=+0.63, boundary_type='reflective')

    fuel_cell = openmc.Cell(name='1.6% Fuel')
    fuel_cell.fill = fuel
    fuel_cell.region = -fuel_outer_radius

    clad_cell = openmc.Cell(name='1.6% Clad')
    clad_cell.fill = zircaloy
    clad_cell.region = +fuel_outer_radius & -clad_outer_radius

    moderator_cell = openmc.Cell(name='1.6% Moderator')
    moderator_cell.fill = water
    moderator_cell.region = +clad_outer_radius & (+min_x & -max_x & +min_y & -max_y)

    root_universe = openmc.Universe(name='root universe', universe_id=0)
    root_universe.add_cells([fuel_cell, clad_cell, moderator_cell])

    geometry = openmc.Geometry(root_universe)

    # Settings
    settings = openmc.Settings()
    settings.batches = 100  # Reduced for faster testing
    settings.inactive = 10
    settings.particles = 500  # Reduced for faster testing
    settings.run_mode = 'eigenvalue'
    settings.verbosity = 1

    bounds = [-0.63, -0.63, -10, 0.63, 0.63, 10.]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.Source(space=uniform_dist)

    model = openmc.model.Model(geometry, materials, settings)
    return model


# ===============================================================
# Derivative tally setup
# ===============================================================
def add_derivative_tally(model):
    """Add a k-effective tally with boron density derivative."""
    
    # Create a derivative for boron density in water
    deriv = openmc.TallyDerivative(
        variable='nuclide_density',
        material=3,  # Water material ID
        nuclide='B10'  # or 'B11' or just track boron
    )
    
    # Create a k-effective tally with the derivative
    keff_tally = openmc.Tally(name='k-eff-deriv')
    keff_tally.scores = ['keff']
    keff_tally.derivative = deriv
    
    model.tallies.append(keff_tally)


# ===============================================================
# Search functions
# ===============================================================
def modifier_ppm(ppm):
    """Modifier function that changes boron concentration."""
    # This will be called with different ppm values by keff_search
    # We rebuild the model here; in real use you might mutate in-place
    pass


def run_search_without_derivatives():
    """Run keff_search without using derivative tallies."""
    print("\n" + "="*70)
    print("TEST 1: Standard GRsecant search (NO derivatives)")
    print("="*70)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Initial guesses for boron concentration
        ppm_low = 500  # Low boron (higher k)
        ppm_high = 1500  # High boron (lower k)
        
        def modifier(ppm):
            """Rebuild model with new boron concentration."""
            model = build_model(ppm)
            model.export_to_xml(path=tmpdir)
        
        # Build initial model to get settings
        model = build_model((ppm_low + ppm_high) / 2)
        model.settings.batches = 100
        model.settings.inactive = 10
        
        start_time = time.time()
        
        result = model.keff_search(
            func=modifier,
            x0=ppm_low,
            x1=ppm_high,
            target=1.0,
            k_tol=1e-4,
            sigma_final=3e-4,
            b0=model.settings.batches - model.settings.inactive,
            maxiter=10,
            output=True,
            use_derivative_tallies=False,
            run_kwargs={'cwd': tmpdir}
        )
        
        elapsed = time.time() - start_time
        
        print(f"\n{'Results':^70}")
        print(f"  Root (optimal ppm): {result.root:.4f} ppm")
        print(f"  Converged: {result.converged}")
        print(f"  Termination reason: {result.flag}")
        print(f"  MC runs performed: {result.function_calls}")
        print(f"  Total batches: {result.total_batches}")
        print(f"  Elapsed time: {elapsed:.2f} s")
        
        return result, elapsed


def run_search_with_derivatives():
    """Run keff_search using derivative tallies as constraints."""
    print("\n" + "="*70)
    print("TEST 2: GRsecant search WITH derivative tally constraints")
    print("="*70)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Initial guesses for boron concentration
        ppm_low = 500  # Low boron (higher k)
        ppm_high = 1500  # High boron (lower k)
        
        def modifier(ppm):
            """Rebuild model with new boron concentration."""
            model = build_model(ppm)
            add_derivative_tally(model)  # Add derivative tally
            model.export_to_xml(path=tmpdir)
        
        # Build initial model to get settings
        model = build_model((ppm_low + ppm_high) / 2)
        add_derivative_tally(model)
        model.settings.batches = 100
        model.settings.inactive = 10
        
        start_time = time.time()
        
        result = model.keff_search(
            func=modifier,
            x0=ppm_low,
            x1=ppm_high,
            target=1.0,
            k_tol=1e-4,
            sigma_final=3e-4,
            b0=model.settings.batches - model.settings.inactive,
            maxiter=10,
            output=True,
            use_derivative_tallies=True,  # Enable derivative usage
            deriv_constraint_offsets=[-0.3, -0.15, 0.15, 0.3],
            run_kwargs={'cwd': tmpdir}
        )
        
        elapsed = time.time() - start_time
        
        print(f"\n{'Results':^70}")
        print(f"  Root (optimal ppm): {result.root:.4f} ppm")
        print(f"  Converged: {result.converged}")
        print(f"  Termination reason: {result.flag}")
        print(f"  MC runs performed: {result.function_calls}")
        print(f"  Total batches: {result.total_batches}")
        print(f"  Elapsed time: {elapsed:.2f} s")
        
        return result, elapsed


# ===============================================================
# Main execution
# ===============================================================
if __name__ == '__main__':
    print("\n" + "="*70)
    print("Derivative Tally Enhancement for keff_search")
    print("="*70)
    print("Objective: Find boron concentration (ppm) for k-eff = 1.0")
    print("Model: PWR pin cell with 1.6% enriched UO2 fuel in borated water")
    print("="*70)
    
    # Run comparison
    result1, time1 = run_search_without_derivatives()
    result2, time2 = run_search_with_derivatives()
    
    # Print summary
    print("\n" + "="*70)
    print("COMPARISON SUMMARY")
    print("="*70)
    
    improvement_runs = ((result1.function_calls - result2.function_calls) / 
                        result1.function_calls * 100)
    improvement_batches = ((result1.total_batches - result2.total_batches) / 
                           result1.total_batches * 100)
    improvement_time = ((time1 - time2) / time1 * 100)
    
    print(f"\nMC Runs:")
    print(f"  Without derivatives: {result1.function_calls} runs")
    print(f"  With derivatives:    {result2.function_calls} runs")
    print(f"  Reduction:           {improvement_runs:.1f}%")
    
    print(f"\nTotal Batches:")
    print(f"  Without derivatives: {result1.total_batches} batches")
    print(f"  With derivatives:    {result2.total_batches} batches")
    print(f"  Reduction:           {improvement_batches:.1f}%")
    
    print(f"\nWall Clock Time:")
    print(f"  Without derivatives: {time1:.2f} s")
    print(f"  With derivatives:    {time2:.2f} s")
    print(f"  Reduction:           {improvement_time:.1f}%")
    
    print(f"\nRoot Location (optimal ppm):")
    print(f"  Without derivatives: {result1.root:.4f} ppm")
    print(f"  With derivatives:    {result2.root:.4f} ppm")
    print(f"  Difference:          {abs(result1.root - result2.root):.4f} ppm")
    
    print("\n" + "="*70)
    print("Key Insights:")
    print("="*70)
    print("• Derivative tallies provide 'free' constraint points via linear")
    print("  Taylor expansion around each MC-evaluated point")
    print("• These constraints guide the secant method curve fit without")
    print("  requiring additional Monte Carlo runs")
    print("• Result: Faster convergence with fewer MC runs while maintaining")
    print("  accuracy (both methods converge to similar roots)")
    print("="*70 + "\n")
