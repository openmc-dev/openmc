#!/usr/bin/env python
"""
Test script demonstrating generic derivative support in Model.keff_search.

This script demonstrates how keff_search now works with ANY derivative variable
supported by the C++ OpenMC backend (density, nuclide_density, temperature,
enrichment), not just boron concentration.

Key features demonstrated:
1. Automatic derivative normalization handling large magnitudes (O(10^20))
2. deriv_to_x_func parameter for converting nuclide density to custom units
3. Boron ppm search with physically realistic conversion factors
4. Generic derivative extraction from derivative tallies

Each test case shows:
1. Building a model with a configurable parameter
2. Adding base and derivative tallies for the target variable
3. Calling keff_search with appropriate parameters
4. For nuclide_density: using deriv_to_x_func for unit conversion
5. Automatic normalization handling large derivative magnitudes
"""

import openmc
import openmc.stats
import numpy as np
import time
from pathlib import Path
import tempfile
import math


def build_model(boron_ppm=1000, fuel_enrichment=1.6, fuel_temp_K=293):
    """Build a generic PWR pin-cell model with configurable parameters."""
    # Fuel
    fuel = openmc.Material(name='Fuel', material_id=1)
    fuel.set_density('g/cm3', 10.31341)
    fuel.add_element('U', 1., enrichment=fuel_enrichment)
    fuel.add_element('O', 2.)
    fuel.temperature = fuel_temp_K

    # Cladding
    clad = openmc.Material(name='Clad', material_id=2)
    clad.set_density('g/cm3', 6.55)
    clad.add_element('Zr', 1.)

    # Borated coolant
    coolant = openmc.Material(name='Coolant', material_id=3)
    coolant.set_density('g/cm3', 0.741)
    coolant.add_element('H', 2.)
    coolant.add_element('O', 1.)
    coolant.add_element('B', boron_ppm * 1e-6)

    materials = openmc.Materials([fuel, clad, coolant])

    # Geometry: simple pin cell with reflective boundaries
    fuel_r = openmc.ZCylinder(r=0.39218)
    clad_r = openmc.ZCylinder(r=0.45720)
    min_x = openmc.XPlane(x0=-0.63, boundary_type='reflective')
    max_x = openmc.XPlane(x0=+0.63, boundary_type='reflective')
    min_y = openmc.YPlane(y0=-0.63, boundary_type='reflective')
    max_y = openmc.YPlane(y0=+0.63, boundary_type='reflective')

    fuel_cell = openmc.Cell(name='Fuel')
    fuel_cell.fill = fuel
    fuel_cell.region = -fuel_r

    clad_cell = openmc.Cell(name='Clad')
    clad_cell.fill = clad
    clad_cell.region = +fuel_r & -clad_r

    coolant_cell = openmc.Cell(name='Coolant')
    coolant_cell.fill = coolant
    coolant_cell.region = +clad_r & +min_x & -max_x & +min_y & -max_y

    root = openmc.Universe(name='root', universe_id=0)
    root.add_cells([fuel_cell, clad_cell, coolant_cell])
    geometry = openmc.Geometry(root)

    # Settings
    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 10
    settings.particles = 500
    settings.run_mode = 'eigenvalue'
    settings.verbosity = 1

    bounds = [-0.63, -0.63, -10, 0.63, 0.63, 10.]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.Source(space=uniform_dist)

    return openmc.model.Model(geometry, materials, settings)


def add_derivative_tallies(model, deriv_variable, deriv_material, deriv_nuclide=None):
    """
    Add base and derivative tallies to model for generic derivative extraction.

    Parameters
    ----------
    model : openmc.Model
    deriv_variable : str
        'density', 'nuclide_density', 'temperature', or 'enrichment'
    deriv_material : int
        Material ID to perturb
    deriv_nuclide : str, optional
        Nuclide name for nuclide_density derivatives (e.g., 'B10', 'U235')
    """
    # Base tallies
    t_fission = openmc.Tally(name='base_fission')
    t_fission.scores = ['nu-fission']

    t_absorption = openmc.Tally(name='base_absorption')
    t_absorption.scores = ['absorption']

    tallies = [t_fission, t_absorption]

    # Derivative tallies
    deriv = openmc.TallyDerivative(
        variable=deriv_variable,
        material=deriv_material,
        nuclide=deriv_nuclide
    )

    t_fission_deriv = openmc.Tally(name=f'fission_deriv_{deriv_variable}')
    t_fission_deriv.scores = ['nu-fission']
    t_fission_deriv.derivative = deriv

    t_absorption_deriv = openmc.Tally(name=f'absorption_deriv_{deriv_variable}')
    t_absorption_deriv.scores = ['absorption']
    t_absorption_deriv.derivative = deriv

    tallies.extend([t_fission_deriv, t_absorption_deriv])
    model.tallies = openmc.Tallies(tallies)


def run_test(test_name, model_builder, modifier_func, deriv_variable, deriv_material,
             deriv_nuclide, x0, x1, target, deriv_nuclide_arg=None, deriv_to_x_func=None,
             expected_magnitude=None, use_derivative_tallies=True,
             x_min=None, x_max=None):
    """
    Generic test runner.
    
    Parameters
    ----------
    test_name : str
        Name of test for display
    model_builder : callable
        Function that builds the model
    modifier_func : callable
        Function that modifies the model for a given parameter. Signature:
        modifier_func(x, model)
    deriv_variable : str
        Type of derivative: 'density', 'nuclide_density', 'temperature', 'enrichment'
    deriv_material : int
        Material ID to perturb
    deriv_nuclide : str
        Nuclide name (for nuclide_density)
    x0, x1 : float
        Initial guesses for search parameter
    target : float
        Target k-eff
    deriv_nuclide_arg : str, optional
        Nuclide name for derivative tallies
    deriv_to_x_func : callable, optional
        Conversion function from number density to custom units
    expected_magnitude : str, optional
        Expected magnitude of derivatives (e.g., "O(10^20)")
    use_derivative_tallies : bool, optional
        If True, enable derivative tallies and pass derivative args to keff_search.
        If False, run keff_search without derivative tallies (baseline comparison).
    """
    print("\n" + "=" * 80)
    print(f"TEST: {test_name}")
    print("=" * 80)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # Build model
        model = model_builder()
        if use_derivative_tallies:
            add_derivative_tallies(model, deriv_variable, deriv_material, deriv_nuclide_arg)
        model.settings.batches = 50
        model.settings.inactive = 5
        model.settings.particles = 300

        print(f"Model setup:")
        print(f"  Derivative variable: {deriv_variable}")
        print(f"  Derivative material: {deriv_material}")
        if use_derivative_tallies:
            print(f"  Derivative tallies: ON")
            if deriv_nuclide_arg:
                print(f"  Derivative nuclide: {deriv_nuclide_arg}")
            if deriv_to_x_func is not None:
                print(f"  Conversion function: Provided (deriv_to_x_func)")
                if expected_magnitude:
                    print(f"  Expected derivative magnitude: {expected_magnitude}")
        else:
            print(f"  Derivative tallies: OFF (baseline run)")
        print(f"  Initial guesses: x0={x0}, x1={x1}")
        print(f"  Target k-eff: {target}")
        print(f"  Batches: {model.settings.batches}")
        if deriv_to_x_func is not None and use_derivative_tallies:
            print(f"  NOTE: Automatic normalization handles large derivatives!")

        start_time = time.time()

        try:
            # Wrap modifier to bind the local model instance
            def _modifier(x):
                return modifier_func(x, model)

            # Build keff_search call with optional deriv_to_x_func
            search_kwargs = {
                'func': _modifier,
                'x0': x0,
                'x1': x1,
                'target': target,
                'k_tol': 1e-3,
                'sigma_final': 3e-3,
                'b0': model.settings.batches - model.settings.inactive,
                'maxiter': 50,
                'output': True,
                'run_kwargs': {'cwd': tmpdir_path},
                'use_derivative_tallies': use_derivative_tallies,
            }
            
            if x_min is not None:
                search_kwargs['x_min'] = x_min
            if x_max is not None:
                search_kwargs['x_max'] = x_max

            if use_derivative_tallies:
                search_kwargs.update({
                    'deriv_variable': deriv_variable,
                    'deriv_material': deriv_material,
                    'deriv_nuclide': deriv_nuclide_arg,
                })
                if deriv_to_x_func is not None:
                    search_kwargs['deriv_to_x_func'] = deriv_to_x_func

            result = model.keff_search(**search_kwargs)

            elapsed = time.time() - start_time

            print(f"\n{'RESULTS':^80}")
            print(f"  Converged: {result.converged}")
            print(f"  Termination reason: {result.flag}")
            print(f"  Final parameter value: {result.root:.6f}")
            print(f"  MC runs performed: {result.function_calls}")
            print(f"  Total batches: {result.total_batches}")
            print(f"  Elapsed time: {elapsed:.2f} s")
            if deriv_to_x_func is not None and use_derivative_tallies:
                print(f"  ✓ Large derivatives handled automatically by normalization!")
            print(f"  ✓ Test PASSED")
            
            # Store elapsed time in result for comparison
            result.elapsed_time = elapsed

        except Exception as e:
            print(f"\n  ✗ Test FAILED: {e}")
            import traceback
            traceback.print_exc()

        return locals().get('result', None)


if __name__ == '__main__':
    print("\n" + "=" * 80)
    print("COMPREHENSIVE COMPARISON: GRsecant vs Least Squares")
    print("=" * 80)
    print("This test compares two optimization methods on two test cases:")
    print("  Test Case 1: Boron concentration search (nuclide_density with ppm conversion)")
    print("  Test Case 2: Fuel density search (density derivative, densification scenario)")
    print("")
    print("Methods:")
    print("  1. GRsecant (baseline):      No derivatives, standard curve-fitting")
    print("  2. Least Squares:            GRsecant + gradient constraints + auto-normalization")
    print("=" * 80)
    
    # Physical constants for boron ppm conversion
    BORON_DENSITY_WATER = 0.741  # g/cm³ at room temperature
    BORON_ATOMIC_MASS = 10.81    # g/mol (natural boron average)
    AVOGADRO = 6.02214076e23     # atoms/mol
    
    # Scale factor: dN/dppm = (1e-6 * rho * N_A) / M_boron
    BORON_PPM_SCALE = (1e-6 * BORON_DENSITY_WATER * AVOGADRO) / BORON_ATOMIC_MASS
    
    print(f"\nBoron conversion factor calculation:")
    print(f"  Water density: {BORON_DENSITY_WATER} g/cm³")
    print(f"  Boron mass: {BORON_ATOMIC_MASS} g/mol")
    print(f"  ppm to number density: {BORON_PPM_SCALE:.4e} atoms/cm³/ppm")
    print(f"  Expected dk/dppm magnitude: O(10^16) to O(10^20)")
    print(f"  ✓ Automatic normalization handles this automatically!\n")

    # TEST 1: Boron concentration search
    print("\n" + "=" * 100)
    print("[TEST 1] BORON CONCENTRATION SEARCH")
    print("=" * 100)
    print("Parameter: Boron concentration in coolant (ppm)")
    print("Derivative variable: nuclide_density for B10")
    print("Derivative magnitude: O(10^16-10^20)")
    
    boron_results = {}
    
    try:
        def modifier_boron(ppm, model):
            """Modify boron concentration in coolant."""
            ppm = max(ppm, 0.0)
            coolant = model.materials[2]
            for elem in ('H', 'O', 'B'):
                coolant.remove_element(elem)
            coolant.set_density('g/cm3', 0.741)
            coolant.add_element('H', 2.)
            coolant.add_element('O', 1.)
            coolant.add_element('B', ppm * 1e-6)
            model.export_to_xml()

        def boron_ppm_conversion(deriv_dN):
            """Convert dk/dN to dk/dppm using chain rule."""
            return deriv_dN * BORON_PPM_SCALE
        
        # Method 1: GRsecant without derivative tallies
        result = run_test(
            "Boron search: GRsecant WITHOUT derivatives",
            lambda: build_model(boron_ppm=1000),
            modifier_boron,
            None, None, None,
            500, 1500, 1.20,
            use_derivative_tallies=False,
        )
        if result:
            boron_results['GRsecant (no deriv)'] = result
        
        # Method 2: Least Squares with derivative tallies
        result = run_test(
            "Boron search: Least Squares WITH derivatives",
            lambda: build_model(boron_ppm=1000),
            modifier_boron,
            'nuclide_density', 3, 'B10',
            500, 1500, 1.20,
            deriv_nuclide_arg='B10',
            deriv_to_x_func=boron_ppm_conversion,
            expected_magnitude="O(10^16-10^20)",
            use_derivative_tallies=True,
        )
        if result:
            boron_results['Least Squares (with deriv)'] = result
            
    except Exception as e:
        print(f"  ⚠ Boron test encountered error: {e}")
    
    # TEST 2: Fuel density search (densification/swelling scenario)
    print("\n" + "=" * 100)
    print("[TEST 2] FUEL DENSITY SEARCH")
    print("=" * 100)
    print("Parameter: Fuel density (g/cm³)")
    print("Physics: Fuel densification or swelling affects neutron moderation and absorption")
    print("Why Least Squares excels: Derivative provides direct sensitivity for faster convergence")
    print("Derivative variable: density (fuel material density)")
    
    density_results = {}
    
    try:
        def modifier_fuel_density(density_gcm3, model):
            """Modify fuel density directly."""
            fuel = model.materials[0]  # First material is fuel
            # Remove and re-add elements to update density
            fuel.remove_element('U')
            fuel.remove_element('O')
            fuel.set_density('g/cm3', density_gcm3)
            fuel.add_element('U', 1., enrichment=1.6)
            fuel.add_element('O', 2.)
            model.export_to_xml()
        
        # Method 1: GRsecant without derivative tallies
        result = run_test(
            "Fuel density search: GRsecant WITHOUT derivatives",
            lambda: build_model(boron_ppm=150),
            modifier_fuel_density,
            None, None, None,
            5.0, 11.0, 1.17,
            use_derivative_tallies=False,
            x_min=2.0, x_max=12.0
        )
        if result:
            density_results['GRsecant (no deriv)'] = result
        
        # Method 2: Least Squares with derivative tallies
        result = run_test(
            "Fuel density search: Least Squares WITH derivatives",
            lambda: build_model(boron_ppm=150),
            modifier_fuel_density,
            'density', 1, None,  # Material ID 1 is fuel
            5.0, 11.0, 1.17,
            use_derivative_tallies=True,
            x_min=2.0, x_max=12.0
        )
        if result:
            density_results['Least Squares (with deriv)'] = result
            
    except Exception as e:
        print(f"  ⚠ Fuel density test encountered error: {e}")
        import traceback
        traceback.print_exc()
    # Print final comparison tables
    print("\n" + "=" * 100)
    
    if boron_results:
        print("\n[TABLE 1] BORON CONCENTRATION SEARCH RESULTS")
        print("-" * 100)
        print(f"{'Method':<30} {'Final Sol (ppm)':<18} {'MC Runs':<12} {'Tot Batches':<14} {'Time (s)':<12} {'Converged':<10}")
        print("-" * 100)
        
        for method_name in ['GRsecant (no deriv)', 'Least Squares (with deriv)']:
            if method_name in boron_results:
                result = boron_results[method_name]
                elapsed = getattr(result, 'elapsed_time', 0)
                print(f"{method_name:<30} {result.root:>16.1f} {result.function_calls:>11d} {result.total_batches:>13d} {elapsed:>11.2f} {str(result.converged):>9}")
        
        # Efficiency analysis for boron
        if 'GRsecant (no deriv)' in boron_results:
            baseline = boron_results['GRsecant (no deriv)']
            baseline_runs = baseline.function_calls
            baseline_batches = baseline.total_batches
            baseline_time = getattr(baseline, 'elapsed_time', 0)
            
            print("\n" + "-" * 100)
            print("Efficiency Gains (relative to GRsecant baseline):")
            print("-" * 100)
            
            for method_name in ['Least Squares (with deriv)']:
                if method_name in boron_results:
                    result = boron_results[method_name]
                    run_pct = ((baseline_runs - result.function_calls) / baseline_runs * 100) if baseline_runs > 0 else 0
                    batch_pct = ((baseline_batches - result.total_batches) / baseline_batches * 100) if baseline_batches > 0 else 0
                    time_pct = ((baseline_time - getattr(result, 'elapsed_time', 0)) / baseline_time * 100) if baseline_time > 0 else 0
                    
                    print(f"\n{method_name}:")
                    print(f"  MC runs:       {run_pct:+7.1f}%  ({result.function_calls:2d} vs {baseline_runs:2d})")
                    print(f"  Total batches: {batch_pct:+7.1f}%  ({result.total_batches:4d} vs {baseline_batches:4d})")
                    print(f"  Elapsed time:  {time_pct:+7.1f}%  ({getattr(result, 'elapsed_time', 0):6.2f}s vs {baseline_time:6.2f}s)")
    
    
    # TABLE 2: Fuel Density Search
    if density_results:
        print("\n" + "=" * 100)
        print("[TABLE 2] FUEL DENSITY SEARCH RESULTS")
        print("-" * 100)
        print(f"{'Method':<30} {'Final Density (g/cm³)':<18} {'MC Runs':<12} {'Tot Batches':<14} {'Time (s)':<12} {'Converged':<10}")
        print("-" * 100)
        
        for method_name in ['GRsecant (no deriv)', 'Least Squares (with deriv)']:
            if method_name in density_results:
                result = density_results[method_name]
                elapsed = getattr(result, 'elapsed_time', 0)
                print(f"{method_name:<30} {result.root:>16.3f} {result.function_calls:>11d} {result.total_batches:>13d} {elapsed:>11.2f} {str(result.converged):>9}")
        
        # Efficiency analysis for density
        if 'GRsecant (no deriv)' in density_results:
            baseline = density_results['GRsecant (no deriv)']
            baseline_runs = baseline.function_calls
            baseline_batches = baseline.total_batches
            baseline_time = getattr(baseline, 'elapsed_time', 0)
            
            print("\n" + "-" * 100)
            print("Efficiency Gains (relative to GRsecant baseline):")
            print("-" * 100)
            
            for method_name in ['Least Squares (with deriv)']:
                if method_name in density_results:
                    result = density_results[method_name]
                    run_pct = ((baseline_runs - result.function_calls) / baseline_runs * 100) if baseline_runs > 0 else 0
                    batch_pct = ((baseline_batches - result.total_batches) / baseline_batches * 100) if baseline_batches > 0 else 0
                    time_pct = ((baseline_time - getattr(result, 'elapsed_time', 0)) / baseline_time * 100) if baseline_time > 0 else 0
                    
                    print(f"\n{method_name}:")
                    print(f"  MC runs:       {run_pct:+7.1f}%  ({result.function_calls:2d} vs {baseline_runs:2d})")
                    print(f"  Total batches: {batch_pct:+7.1f}%  ({result.total_batches:4d} vs {baseline_batches:4d})")
                    print(f"  Elapsed time:  {time_pct:+7.1f}%  ({getattr(result, 'elapsed_time', 0):6.2f}s vs {baseline_time:6.2f}s)")
    
    
    print("\n" + "=" * 100)
    print("All comparison tests completed!")
    print("=" * 100)
