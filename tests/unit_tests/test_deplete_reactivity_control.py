""" Tests for ReactivityController class """

from pathlib import Path

import pytest
import numpy as np

import openmc
import openmc.lib
from openmc.deplete import CoupledOperator
from openmc.deplete import ReactivityController

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"

@pytest.fixture
def model():
    f = openmc.Material(name="fuel")
    f.add_element("U", 1, percent_type="ao", enrichment=4.25)
    f.add_element("O", 2)
    f.set_density("g/cc", 10.4)
    f.temperature = 293.15

    w = openmc.Material(name="water")
    w.add_element("O", 1)
    w.add_element("H", 2)
    w.set_density("g/cc", 1.0)
    w.temperature = 293.15
    w.depletable = True

    h = openmc.Material(name='helium')
    h.add_element('He', 1)
    h.set_density('g/cm3', 0.001598)

    radii = [0.42, 0.45]
    height = 0.5

    f.volume = np.pi * radii[0] ** 2 * height
    w.volume = np.pi * (radii[1]**2 - radii[0]**2) * height/2

    materials = openmc.Materials([f, w, h])

    surf_interface = openmc.ZPlane(z0=0)
    surf_top = openmc.ZPlane(z0=height/2)
    surf_bot = openmc.ZPlane(z0=-height/2)
    surf_in = openmc.Sphere(r=radii[0])
    surf_out = openmc.Sphere(r=radii[1], boundary_type='vacuum')

    cell_water = openmc.Cell(fill=w, region=-surf_interface)
    cell_helium = openmc.Cell(fill=h, region=+surf_interface)
    universe = openmc.Universe(cells=(cell_water, cell_helium))
    cell_fuel = openmc.Cell(name='fuel_cell', fill=f,
                            region=-surf_in & -surf_top & +surf_bot)
    cell_universe = openmc.Cell(name='universe_cell',fill=universe,
                            region=+surf_in & -surf_out & -surf_top & +surf_bot)
    geometry = openmc.Geometry([cell_fuel, cell_universe])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)

@pytest.fixture
def operator(model):
    return CoupledOperator(model, CHAIN_PATH)

@pytest.fixture
def integrator(operator):
    return openmc.deplete.PredictorIntegrator(
            operator, [1,1], 0.0, timestep_units = 'd')

def test_reactivity_controller_init(operator):
    """Test ReactivityController initialization"""
    def dummy_function(x):
        return x
    
    # Valid initialization
    controller = ReactivityController(
        operator=operator,
        function=dummy_function,
        x0=0.0,
        x1=1.0,
        bracket=[0.0, 2.0]
    )
    
    assert controller.operator == operator
    assert controller.function == dummy_function
    assert controller.x0 == 0.0
    assert controller.x1 == 1.0
    assert controller.kwargs['x_min'] == 0.0
    assert controller.kwargs['x_max'] == 2.0
    
    # Test invalid bracket with wrong length
    with pytest.raises(ValueError, match="bracket must have exactly 2 elements"):
        ReactivityController(operator, dummy_function, 0.0, 1.0, [0.0])
    
    # Test invalid bracket with wrong order
    with pytest.raises(ValueError, match="bracket\\[0\\] must be < bracket\\[1\\]"):
        ReactivityController(operator, dummy_function, 0.0, 1.0, [2.0, 1.0])

def test_reactivity_controller_call(operator):
    """Test the __call__ method acts as a proxy to the function"""
    call_log = []
    
    def test_function(x):
        call_log.append(x)
        return x * 2
    
    controller = ReactivityController(
        operator=operator,
        function=test_function,
        x0=0.0,
        x1=1.0,
        bracket=[0.0, 2.0]
    )
    
    # Test that calling the controller calls the underlying function
    result = controller(5.0)
    assert result == 10.0
    assert call_log == [5.0]
    
    # Test with multiple arguments
    def multi_arg_function(x, y, z=3):
        return x + y + z
    
    controller2 = ReactivityController(
        operator=operator,
        function=multi_arg_function,
        x0=0.0,
        x1=1.0,
        bracket=[0.0, 2.0]
    )
    
    result = controller2(1, 2, z=4)
    assert result == 7

def translate_cell(position):
    """Helper function to translate a cell"""
    cell = [c for c in openmc.lib.cells.values() if c.name == 'universe_cell'][0]
    openmc.lib.cells[cell.id].translation = [0, 0, position]
    return position


def rotate_cell(angle):
    """Helper function to rotate a cell"""
    cell = [c for c in openmc.lib.cells.values() if c.name == 'universe_cell'][0]
    openmc.lib.cells[cell.id].rotation = [0, 0, angle]
    return angle


def adjust_fuel_density(density_factor):
    """Helper function to adjust fuel density"""
    fuel = [m for m in openmc.lib.materials.values() if m.name == 'fuel'][0]
    nuclides = openmc.lib.materials[fuel.id].nuclides
    current_densities = openmc.lib.materials[fuel.id].densities
    new_densities = [d * density_factor for d in current_densities]
    openmc.lib.materials[fuel.id].set_densities(nuclides, new_densities)
    return density_factor

@pytest.mark.parametrize("function, x0, x1, bracket, test_value", [
    (translate_cell, -1.0, 1.0, [-5.0, 5.0], 0.5),
    (rotate_cell, -45.0, 45.0, [-90.0, 90.0], 10.0),
    (adjust_fuel_density, 0.8, 1.2, [0.5, 1.5], 1.0)
])
def test_reactivity_controller_with_openmc_functions(
    run_in_tmpdir, model, operator, function, x0, x1, bracket, test_value
):
    """Test ReactivityController with actual OpenMC geometry modification functions"""
    
    controller = ReactivityController(
        operator=operator,
        function=function,
        x0=x0,
        x1=x1,
        bracket=bracket,
        kwargs={'output': False, 'k_tol': 0.1}
    )
    
    # Export model and initialize OpenMC library
    model.export_to_xml()
    openmc.lib.init()
    
    try:
        # Test that the controller can be called to modify the geometry
        result = controller(test_value)
        assert result == test_value
        
        # Verify the function was actually called by checking it doesn't raise
        # (actual verification would require checking geometry state)
        
    finally:
        openmc.lib.finalize()

def test_controller_kwargs_handling(operator):
    """Test that additional kwargs are properly stored and accessible"""
    def dummy_func(x):
        return x
    
    custom_kwargs = {
        'target': 1.0,
        'k_tol': 1e-4,
        'output': True
    }
    
    controller = ReactivityController(
        operator=operator,
        function=dummy_func,
        x0=0.0,
        x1=1.0,
        bracket=[0.0, 2.0],
        kwargs=custom_kwargs
    )
    
    # Check that bracket limits were added to kwargs
    assert controller.kwargs['x_min'] == 0.0
    assert controller.kwargs['x_max'] == 2.0
    
    # Check that custom kwargs were preserved
    assert controller.kwargs['target'] == 1.0
    assert controller.kwargs['k_tol'] == 1e-4
    assert controller.kwargs['output'] == True

def test_integrator_add_reactivity_control(run_in_tmpdir, model, operator, integrator):
    """Test adding reactivity control to integrator using the add_reactivity_control method"""
    
    def adjust_density(factor):
        # This would modify material density in practice
        return factor
    
    # Test that add_reactivity_control method exists and works
    integrator.add_reactivity_control(
        function=adjust_density,
        x0=0.8,
        x1=1.2,
        bracket=[0.5, 1.5],
        kwargs={'k_tol': 0.1, 'output': False}
    )
    
    assert integrator.reactivity_control is not None
    assert isinstance(integrator.reactivity_control, ReactivityController)
    assert integrator.reactivity_control.x0 == 0.8
    assert integrator.reactivity_control.x1 == 1.2
    assert integrator.reactivity_control.kwargs['x_min'] == 0.5
    assert integrator.reactivity_control.kwargs['x_max'] == 1.5   

