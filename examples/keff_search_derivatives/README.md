# keff_search with Derivative Tallies

This example demonstrates the derivative-accelerated k-effective search capability in OpenMC, which uses gradient information from derivative tallies to significantly speed up convergence during criticality searches.

## Overview

The `Model.keff_search()` method can leverage derivative tallies to compute sensitivities (dk/dx) with respect to material properties, enabling faster and more robust convergence compared to traditional derivative-free methods. This feature is inspired by the methodology developed by Sterling Harper in "Calculating Reaction Rate Derivatives in Monte Carlo Neutron Transport" (https://dspace.mit.edu/bitstream/handle/1721.1/106690/969775837-MIT.pdf).

This example compares:

1. **GRsecant (baseline)**: Standard gradient-free search using only k-effective values
2. **Least Squares with Derivatives**: Enhanced search using both k-effective values and derivative constraints

## Key Features

- **Automatic derivative tally setup**: No manual tally configuration required
- **Automatic derivative normalization**: Handles derivatives with very large magnitudes (e.g., O(10²⁰) for ppm-scale derivatives)
- **Generic derivative support**: Works with supported derivative variables:
  - `nuclide_density`: Perturbations to specific nuclide concentrations ✓ **Fully supported**
  - `density`: Material mass density changes ✓ **Fully supported**
  - `temperature`: Doppler temperature effects ⚠️ **Limited support** (see below)

### Temperature Derivative Limitations

Temperature derivatives in OpenMC have **significant limitations** that make them impractical for most k-eff search applications:

1. **Requires Windowed Multipole (WMP) data**: Temperature derivatives are only computed when using multipole cross section representation. Standard tabulated cross sections (HDF5, ACE) at discrete temperatures do not support analytical temperature derivatives.

2. **Limited energy range**: Derivatives are only valid within the resolved resonance energy range where the multipole approximation is applicable (~1 eV to ~10 keV for most materials). Outside this range, derivatives return zero.

3. **Not implemented for interpolated cross sections**: If using OpenMC's temperature interpolation feature (`settings.temperature_method = 'interpolation'`), temperature derivatives are not computed from the interpolation - only from multipole data.

4. **Sparse multipole data availability**: Very few nuclides in standard nuclear data libraries include WMP data. Most reactor applications would have insufficient coverage.

**Recommendation**: For k-eff searches involving temperature changes, use the standard derivative-free search (set `use_derivative_tallies=False`) rather than attempting to use temperature derivatives. The limitations above make temperature derivatives unreliable for reactor-scale criticality searches.

## Files

- `test_tally_deriv_keff_search.py`: Comprehensive comparison script demonstrating two search problems:
  1. **Boron concentration search**: Finding critical boron concentration in PWR coolant (nuclide_density derivative)
  2. **Fuel density search**: Finding critical fuel density for densification scenarios (density derivative)

## Quick Start

### Basic Usage

```bash
python test_tally_deriv_keff_search.py
```

This will run both test cases and display comparison results showing:
- Final converged parameter values
- Number of Monte Carlo runs required
- Total batches executed
- Elapsed time
- Efficiency gains relative to baseline

### Expected Output

The script displays:
1. Physical constants and conversion factors
2. Progress for each iteration (parameter value, k-effective, derivative)
3. Final results table comparing GRsecant vs Least Squares
4. Efficiency analysis showing speedup and resource savings

## Test Cases

### Test Case 1: Boron Concentration Search

**Problem**: Find the boron concentration (in ppm) in PWR coolant to achieve a target k-effective of 1.20

**Physics**: 
- Boron-10 is a strong thermal neutron absorber
- Small changes in concentration significantly affect reactivity

**Derivative Type**: `nuclide_density` for B10

**Key Challenge**: Derivatives are O(10¹⁶-10²⁰) due to unit conversion from atoms/cm³ to ppm. Automatic normalization handles this transparently.

### Test Case 2: Fuel Density Search

**Problem**: Find the fuel density (g/cm³) to achieve a target k-effective of 1.17

**Physics**:
- Fuel densification (aging) or swelling affects neutron moderation
- Density changes affect both fission rates and neutron leakage

**Derivative Type**: `density` for fuel material

**Key Advantage**: Derivatives are O(1), making gradient information highly effective without unit conversion complexity.

## Understanding the Results

### When Derivatives Help Most

1. **Non-linear relationships**: When parameter-to-keff mapping is complex
2. **Large search ranges**: When initial guesses are far from solution
3. **High-precision requirements**: When tight convergence tolerances are needed
4. **Expensive evaluations**: When each MC run is computationally costly

### Convergence Indicators

The script prints iteration details:
```
Iteration 1: batches=45, x=500, keff=1.15234 +/- 0.00234, dk/dx=2.3e+16
```

- `x`: Current parameter value
- `keff`: Computed k-effective with uncertainty
- `dk/dx`: Derivative (sensitivity) extracted from tallies

## Implementation Details

### Automatic Derivative Tally Setup

When `use_derivative_tallies=True`, the `keff_search` method automatically:
1. Adds base tallies (nu-fission, absorption)
2. Adds derivative tallies for the specified variable
3. Extracts derivatives from statepoint files
4. Normalizes derivatives for numerical stability
5. Incorporates derivative constraints into least-squares fitting

**No manual tally configuration required!**

### Derivative Normalization

Large derivatives (e.g., dk/dppm ~ 10²⁰) are automatically normalized using the geometric mean of recent derivative magnitudes:

```
deriv_scale = exp(mean(log(|dk/dx|)))
```

This ensures numerical stability without requiring manual scaling by the user.

### Unit Conversion

For `nuclide_density` derivatives, OpenMC computes dk/dN (where N is number density in atoms/cm³). If your search parameter is in different units (e.g., ppm), provide a conversion function:

```python
deriv_to_x_func=lambda deriv: deriv * conversion_factor
```

For other derivative types (`density`, `temperature`), no conversion is typically needed.

## API Summary

### Minimal Code Pattern

```python
import openmc
import openmc.model

# 1. Build a model with configurable parameters
def build_model(boron_ppm=1000):
    # Define materials
    fuel = openmc.Material(material_id=1)
    fuel.set_density('g/cm3', 10.31341)
    fuel.add_element('U', 1., enrichment=1.6)
    fuel.add_element('O', 2.)
    
    coolant = openmc.Material(material_id=3)
    coolant.set_density('g/cm3', 0.741)
    coolant.add_element('H', 2.)
    coolant.add_element('O', 1.)
    coolant.add_element('B', boron_ppm * 1e-6)
    
    # Define geometry and settings...
    # (see test_tally_deriv_keff_search.py for complete example)
    
    return openmc.model.Model(geometry, materials, settings)

# 2. Define parameter modifier function
def set_boron_ppm(ppm, model):
    """Modify boron concentration in coolant."""
    coolant = model.materials[2]
    coolant.remove_element('H')
    coolant.remove_element('O')
    coolant.remove_element('B')
    coolant.set_density('g/cm3', 0.741)
    coolant.add_element('H', 2.)
    coolant.add_element('O', 1.)
    coolant.add_element('B', ppm * 1e-6)
    model.export_to_xml()

# 3. Run search with derivative tallies (added automatically!)
model = build_model(boron_ppm=1000)

result = model.keff_search(
    func=lambda x: set_boron_ppm(x, model),
    x0=500.0,
    x1=1500.0,
    target=1.20,
    k_tol=1e-3,
    sigma_final=3e-3,
    maxiter=10,
    use_derivative_tallies=True,
    deriv_variable='nuclide_density',
    deriv_material=3,
    deriv_nuclide='B10'
)
```

### Key Parameters

- `use_derivative_tallies` (bool): Enable derivative-accelerated search
- `deriv_variable` (str): `'density'`, `'nuclide_density'`, or `'temperature'`
- `deriv_material` (int): Material ID to perturb
- `deriv_nuclide` (str): Nuclide name (required for `nuclide_density`)
- `deriv_to_x_func` (callable): Unit conversion function (optional)


