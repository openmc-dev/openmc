# keff_search with Derivative Tallies

This example demonstrates the derivative-accelerated k-effective search capability in OpenMC, which uses gradient information from derivative tallies to significantly speed up convergence during criticality searches.

## Overview

The `Model.keff_search()` method can leverage derivative tallies to compute sensitivities (dk/dx) with respect to material properties, enabling faster and more robust convergence compared to traditional derivative-free methods. This example compares:

1. **GRsecant (baseline)**: Standard gradient-free search using only k-effective values
2. **Least Squares with Derivatives**: Enhanced search using both k-effective values and derivative constraints

## Key Features

- **Automatic derivative tally setup**: No manual tally configuration required
- **Automatic derivative normalization**: Handles derivatives with very large magnitudes (e.g., O(10²⁰) for ppm-scale derivatives)
- **Generic derivative support**: Works with any derivative variable supported by OpenMC:
  - `nuclide_density`: Perturbations to specific nuclide concentrations
  - `density`: Material mass density changes
  - `temperature`: Doppler temperature effects

## Files

- `test_generic_keff_search.py`: Comprehensive comparison script demonstrating two search problems:
  1. **Boron concentration search**: Finding critical boron concentration in PWR coolant (nuclide_density derivative)
  2. **Fuel density search**: Finding critical fuel density for densification scenarios (density derivative)

## Requirements

- OpenMC with derivative tally support (C++ backend must be compiled with derivative capability)
- Nuclear cross section data (NNDC HDF5 library recommended)
- Python packages: `numpy`, `scipy`, `openmc`

## Quick Start

### Basic Usage

```bash
python test_generic_keff_search.py
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
- Typical PWR operational range: 0-2000 ppm

**Derivative Type**: `nuclide_density` for B10

**Key Challenge**: Derivatives are O(10¹⁶-10²⁰) due to unit conversion from atoms/cm³ to ppm. Automatic normalization handles this transparently.

**Usage Example**:
```python
model = build_model(boron_ppm=1000)

def set_boron_ppm(ppm):
    coolant = model.materials[2]
    coolant.remove_element('B')
    coolant.add_element('B', ppm * 1e-6)

# Conversion factor: ppm to atoms/cm³
scale = 1e-6 * 0.741 * 6.022e23 / 10.81  # ~4.1e16

result = model.keff_search(
    set_boron_ppm, 500, 1500, target=1.20,
    use_derivative_tallies=True,
    deriv_variable='nuclide_density',
    deriv_material=3,
    deriv_nuclide='B10',
    deriv_to_x_func=lambda d: d * scale
)
```

### Test Case 2: Fuel Density Search

**Problem**: Find the fuel density (g/cm³) to achieve a target k-effective of 1.17

**Physics**:
- Fuel densification (aging) or swelling affects neutron moderation
- Density changes affect both fission rates and neutron leakage
- Typical UO₂ density: 10-11 g/cm³

**Derivative Type**: `density` for fuel material

**Key Advantage**: Derivatives are O(1), making gradient information highly effective without unit conversion complexity.

**Usage Example**:
```python
model = build_model()

def set_fuel_density(density_gcm3):
    fuel = model.materials[0]
    fuel.set_density('g/cm3', density_gcm3)

result = model.keff_search(
    set_fuel_density, 5.0, 11.0, target=1.17,
    use_derivative_tallies=True,
    deriv_variable='density',
    deriv_material=1,
    x_min=2.0, x_max=12.0
)
```

## Understanding the Results

### Typical Performance Gains

Using derivative tallies typically provides:
- **30-50% fewer MC runs**: Fewer iterations to convergence
- **20-40% fewer total batches**: More efficient batch allocation
- **Faster wall-clock time**: Reduced computational cost

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
# 1. Create model
model = openmc.examples.pwr_pin_cell()

# 2. Define parameter modifier
def modify_parameter(x):
    # Modify model based on parameter x
    material.set_property(x)

# 3. Run search (tallies added automatically!)
result = model.keff_search(
    modify_parameter,
    x0, x1,                          # Initial guesses
    target=1.0,                      # Target k-effective
    use_derivative_tallies=True,     # Enable derivatives
    deriv_variable='density',        # Type of derivative
    deriv_material=1,                # Material ID
    deriv_nuclide='B10'              # (if nuclide_density)
)
```

### Key Parameters

- `use_derivative_tallies` (bool): Enable derivative-accelerated search
- `deriv_variable` (str): `'density'`, `'nuclide_density'`, or `'temperature'`
- `deriv_material` (int): Material ID to perturb
- `deriv_nuclide` (str): Nuclide name (required for `nuclide_density`)
- `deriv_to_x_func` (callable): Unit conversion function (optional)

## Customization

### Adjusting Search Parameters

You can tune the search behavior:

```python
result = model.keff_search(
    func, x0, x1, target=1.0,
    k_tol=1e-3,              # Convergence tolerance on k-effective
    sigma_final=3e-3,        # Maximum accepted uncertainty
    maxiter=50,              # Maximum iterations
    b0=90,                   # Initial number of active batches
    b_min=20,                # Minimum batches per iteration
    b_max=200,               # Maximum batches per iteration
    memory=4,                # Points used in curve fitting
    output=True,             # Print iteration details
    use_derivative_tallies=True,
    deriv_variable='nuclide_density',
    deriv_material=3,
    deriv_nuclide='B10'
)
```

### Custom Model Setup

Replace `build_model()` with your own model constructor:

```python
def my_custom_model():
    # Define materials
    fuel = openmc.Material(...)
    
    # Define geometry
    geometry = openmc.Geometry(...)
    
    # Define settings
    settings = openmc.Settings(...)
    
    return openmc.Model(geometry, materials, settings)
```

## Troubleshooting

### "No cross_sections.xml file found"

Set the `OPENMC_CROSS_SECTIONS` environment variable:
```bash
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml
```

Or download NNDC data:
```bash
wget -q -O - https://anl.box.com/shared/static/teaup95cqv8s9nn56hfn7ku8mmelr95p.xz | tar -C $HOME -xJ
export OPENMC_CROSS_SECTIONS=$HOME/nndc_hdf5/cross_sections.xml
```

### Derivatives Not Found

If derivatives are missing from statepoint files:
1. Ensure OpenMC was compiled with derivative support
2. Check that `deriv_variable`, `deriv_material`, and `deriv_nuclide` match your model
3. Verify the C++ backend supports the requested derivative type

### Slow Convergence

If convergence is slower than expected:
1. Check that initial guesses (`x0`, `x1`) bracket the solution
2. Increase `b0` (initial batches) for better statistics
3. Verify derivative magnitudes are reasonable (check iteration output)
4. Try adjusting `deriv_weight` parameter (default is 1.0)

## References

- Price, D., & Roskoff, N. (2023). "GRsecant: A root-finding algorithm with uncertainty quantification for Monte Carlo simulations." *Progress in Nuclear Energy*, 104731.
- OpenMC Documentation: https://docs.openmc.org/
- OpenMC Derivative Tallies: See `include/openmc/tallies/derivative.h` in the OpenMC source

## Contributing

To extend this example:
1. Add new test cases for different derivative types (e.g., temperature)
2. Demonstrate multi-parameter searches (currently single-parameter only)
3. Add visualization of convergence behavior
4. Include more complex geometries (assemblies, full core)

## License

This example is part of OpenMC and is distributed under the MIT License. See the main OpenMC repository for details.
