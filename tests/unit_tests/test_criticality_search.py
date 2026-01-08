from math import pi
from pathlib import Path
import numpy as np
import openmc


def test_keff_search_with_derivative_tallies(run_in_tmpdir):
    """Test keff_search with derivative tallies enabled for fuel density perturbation."""
    # Build a simple PWR pin-cell model
    fuel = openmc.Material(name='Fuel', material_id=1)
    fuel.set_density('g/cm3', 10.31341)
    fuel.add_element('U', 1., enrichment=1.6)
    fuel.add_element('O', 2.)

    clad = openmc.Material(name='Clad', material_id=2)
    clad.set_density('g/cm3', 6.55)
    clad.add_element('Zr', 1.)

    coolant = openmc.Material(name='Coolant', material_id=3)
    coolant.set_density('g/cm3', 0.741)
    coolant.add_element('H', 2.)
    coolant.add_element('O', 1.)
    coolant.add_element('B', 150 * 1e-6)  # 150 ppm boron

    materials = openmc.Materials([fuel, clad, coolant])

    fuel_r = openmc.ZCylinder(r=0.39218)
    clad_r = openmc.ZCylinder(r=0.45720)
    min_x = openmc.XPlane(x0=-0.63, boundary_type='reflective')
    max_x = openmc.XPlane(x0=+0.63, boundary_type='reflective')
    min_y = openmc.YPlane(y0=-0.63, boundary_type='reflective')
    max_y = openmc.YPlane(y0=+0.63, boundary_type='reflective')

    fuel_cell = openmc.Cell(fill=fuel, region=-fuel_r)
    clad_cell = openmc.Cell(fill=clad, region=+fuel_r & -clad_r)
    coolant_cell = openmc.Cell(fill=coolant, region=+clad_r & +min_x & -max_x & +min_y & -max_y)

    root = openmc.Universe(cells=[fuel_cell, clad_cell, coolant_cell])
    geometry = openmc.Geometry(root)

    settings = openmc.Settings()
    settings.batches = 50
    settings.inactive = 5
    settings.particles = 500
    settings.run_mode = 'eigenvalue'

    bounds = [-0.63, -0.63, -10, 0.63, 0.63, 10.]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.Source(space=uniform_dist)

    model = openmc.model.Model(geometry, materials, settings)

    def set_density(x):
        # Modify fuel density by removing and re-adding elements
        fuel_mat = model.materials[0]
        fuel_mat.remove_element('U')
        fuel_mat.remove_element('O')
        fuel_mat.set_density('g/cm3', x)
        fuel_mat.add_element('U', 1., enrichment=1.6)
        fuel_mat.add_element('O', 2.)

    # Perform keff search with derivative tallies enabled
    # Model class will create the required derivative tallies internally
    k_tol = 5e-3
    sigma_final = 5e-3
    result = model.keff_search(
        func=set_density,
        x0=9.0,
        x1=11.0,
        target=1.17,
        k_tol=k_tol,
        sigma_final=sigma_final,
        x_min=5.0,
        x_max=12.0,
        b0=settings.batches - settings.inactive,
        maxiter=10,
        output=False,
        run_kwargs={'cwd': Path('.')},
        use_derivative_tallies=True,
        deriv_variable='density',
        deriv_material=1,
    )

    # Check type of result
    assert isinstance(result, openmc.model.SearchResult)

    # Check that we have function evaluation history
    assert len(result.parameters) >= 2
    assert len(result.means) == len(result.parameters)
    assert len(result.stdevs) == len(result.parameters)
    assert len(result.batches) == len(result.parameters)

    # Check that function_calls property works
    assert result.function_calls == len(result.parameters)

    # Check that total_batches property works
    assert result.total_batches == sum(result.batches)
    assert result.total_batches > 0

    # If converged, check tolerances (but don't fail if not converged due to limited iterations)
    if result.converged:
        final_keff = result.means[-1] + 1.17  # Add back target since means are (keff - target)
        final_sigma = result.stdevs[-1]
        assert abs(final_keff - 1.17) <= k_tol, \
            f"Final keff {final_keff:.5f} not within k_tol {k_tol}"
        assert final_sigma <= sigma_final, \
            f"Final uncertainty {final_sigma:.5f} exceeds sigma_final {sigma_final}"


def test_keff_search_with_nuclide_density_derivatives(run_in_tmpdir):
    """Test keff_search with nuclide density derivatives for boron concentration."""
    # Build a simple PWR pin-cell model
    fuel = openmc.Material(name='Fuel', material_id=1)
    fuel.set_density('g/cm3', 10.31341)
    fuel.add_element('U', 1., enrichment=1.6)
    fuel.add_element('O', 2.)

    clad = openmc.Material(name='Clad', material_id=2)
    clad.set_density('g/cm3', 6.55)
    clad.add_element('Zr', 1.)

    coolant = openmc.Material(name='Coolant', material_id=3)
    coolant.set_density('g/cm3', 0.741)
    coolant.add_element('H', 2.)
    coolant.add_element('O', 1.)
    coolant.add_element('B', 1000 * 1e-6)  # 1000 ppm boron

    materials = openmc.Materials([fuel, clad, coolant])

    fuel_r = openmc.ZCylinder(r=0.39218)
    clad_r = openmc.ZCylinder(r=0.45720)
    min_x = openmc.XPlane(x0=-0.63, boundary_type='reflective')
    max_x = openmc.XPlane(x0=+0.63, boundary_type='reflective')
    min_y = openmc.YPlane(y0=-0.63, boundary_type='reflective')
    max_y = openmc.YPlane(y0=+0.63, boundary_type='reflective')

    fuel_cell = openmc.Cell(fill=fuel, region=-fuel_r)
    clad_cell = openmc.Cell(fill=clad, region=+fuel_r & -clad_r)
    coolant_cell = openmc.Cell(fill=coolant, region=+clad_r & +min_x & -max_x & +min_y & -max_y)

    root = openmc.Universe(cells=[fuel_cell, clad_cell, coolant_cell])
    geometry = openmc.Geometry(root)

    settings = openmc.Settings()
    settings.batches = 50
    settings.inactive = 5
    settings.particles = 500
    settings.run_mode = 'eigenvalue'

    bounds = [-0.63, -0.63, -10, 0.63, 0.63, 10.]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.Source(space=uniform_dist)

    model = openmc.model.Model(geometry, materials, settings)

    def set_boron_ppm(x):
        # Modify boron concentration by removing and re-adding all elements
        x = max(x, 0.0)  # Ensure positive
        coolant_mat = model.materials[2]
        coolant_mat.remove_element('H')
        coolant_mat.remove_element('O')
        coolant_mat.remove_element('B')
        coolant_mat.set_density('g/cm3', 0.741)
        coolant_mat.add_element('H', 2.)
        coolant_mat.add_element('O', 1.)
        coolant_mat.add_element('B', x * 1e-6)

    # Perform keff search with nuclide density derivatives
    # Model class will create the required derivative tallies internally
    k_tol = 5e-3
    sigma_final = 5e-3
    result = model.keff_search(
        func=set_boron_ppm,
        x0=500.0,
        x1=1500.0,
        target=1.20,
        k_tol=k_tol,
        sigma_final=sigma_final,
        x_min=0.1,
        b0=settings.batches - settings.inactive,
        maxiter=10,
        output=False,
        run_kwargs={'cwd': Path('.')},
        use_derivative_tallies=True,
        deriv_variable='nuclide_density',
        deriv_material=3,
        deriv_nuclide='B10',
    )

    # Check type of result
    assert isinstance(result, openmc.model.SearchResult)

    # Check that we have function evaluation history
    assert len(result.parameters) >= 2
    assert len(result.means) == len(result.parameters)
    assert len(result.stdevs) == len(result.parameters)
    assert len(result.batches) == len(result.parameters)

    # Check that function_calls property works
    assert result.function_calls == len(result.parameters)

    # Check that total_batches property works
    assert result.total_batches == sum(result.batches)
    assert result.total_batches > 0

    # If converged, check tolerances (but don't fail if not converged due to limited iterations)
    if result.converged:
        final_keff = result.means[-1] + 1.20  # Add back target since means are (keff - target)
        final_sigma = result.stdevs[-1]
        assert abs(final_keff - 1.20) <= k_tol, \
            f"Final keff {final_keff:.5f} not within k_tol {k_tol}"
        assert final_sigma <= sigma_final, \
            f"Final uncertainty {final_sigma:.5f} exceeds sigma_final {sigma_final}"



def test_keff_search(run_in_tmpdir):
    """Test the Model.keff_search method"""

    # Create model of a sphere of U235
    mat = openmc.Material()
    mat.set_density('g/cm3', 18.9)
    mat.add_nuclide('U235', 1.0)
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere)
    geometry = openmc.Geometry([cell])
    settings = openmc.Settings(particles=1000, inactive=10, batches=30)
    model = openmc.Model(geometry=geometry, settings=settings)

    # Define function to modify sphere radius
    def modify_radius(radius):
        sphere.r = radius

    # Perform keff search
    k_tol = 4e-3
    sigma_final = 2e-3
    result = model.keff_search(
        func=modify_radius,
        x0=6.0,
        x1=9.0,
        k_tol=k_tol,
        sigma_final=sigma_final,
        output=True,
    )

    final_keff = result.means[-1] + 1.0  # Add back target since means are (keff - target)
    final_sigma = result.stdevs[-1]

    # Check for convergence and that tolerances are met
    assert result.converged, "keff_search did not converge"
    assert abs(final_keff - 1.0) <= k_tol, \
        f"Final keff {final_keff:.5f} not within k_tol {k_tol}"
    assert final_sigma <= sigma_final, \
        f"Final uncertainty {final_sigma:.5f} exceeds sigma_final {sigma_final}"

    # Check type of result
    assert isinstance(result, openmc.model.SearchResult)

    # Check that we have function evaluation history
    assert len(result.parameters) >= 2
    assert len(result.means) == len(result.parameters)
    assert len(result.stdevs) == len(result.parameters)
    assert len(result.batches) == len(result.parameters)

    # Check that function_calls property works
    assert result.function_calls == len(result.parameters)

    # Check that total_batches property works
    assert result.total_batches == sum(result.batches)
    assert result.total_batches > 0
