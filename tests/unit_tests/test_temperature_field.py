"""Test the temperature field Python object."""

import openmc
import pytest
import openmc.lib
from openmc.examples import pwr_pin_cell


@pytest.fixture(scope="module")
def mesh():
    """2x2x2 regular mesh"""
    dim = 2
    lower_left = (0., 0., 0.)
    upper_right = (5.0, 5.0, 5.0)
    mesh = openmc.RegularMesh()
    mesh.lower_left = lower_left
    mesh.upper_right = upper_right
    mesh.dimension = (dim, dim, dim)
    return mesh


def test_xml_serialization(mesh, run_in_tmpdir):
    """Test XMl serialization for a temperature field declaration."""
    # Define values
    values = [0., 1., 2., 3., 4., 5., 6., 7.]

    # Create field
    temperature_field = openmc.TemperatureField(mesh, values)

    # Export settings
    settings = openmc.Settings()
    settings.temperature_field = temperature_field
    settings.export_to_xml()

    # Read settings from xml file
    read_settings = openmc.Settings.from_xml()

    # Check consistency
    assert read_settings.temperature_field == temperature_field


def test_invalid_mesh_type():
    """Check one non-implemented mesh type."""
    mesh = openmc.SphericalMesh(r_grid=[0., 1.0])
    values = []

    with pytest.raises(TypeError):
        openmc.TemperatureField(mesh, values)


def test_inconsistent_values_and_cells_sizes(mesh):
    """Check inconsistency in number of values / mesh cells."""
    values = [0.0] * 7

    with pytest.raises(ValueError):
        openmc.TemperatureField(mesh, values)


def test_invalid_values(mesh):
    """Check negative values."""
    values = [-1.0] * 8

    with pytest.raises(ValueError):
        openmc.TemperatureField(mesh, values)


def test_non_numeric_values(mesh):
    """Check that field values are numeric."""
    values = ['not a temperature'] * 8

    with pytest.raises(TypeError):
        openmc.TemperatureField(mesh, values)


def test_lib_temperature_field():
    """Check the Python bindings for the temperature-field C API."""
    model = pwr_pin_cell()
    mesh = openmc.RegularMesh()
    mesh.dimension = (1, 1, 1)
    mesh.lower_left = (-1.0, -1.0, -1.0)
    mesh.upper_right = (1.0, 1.0, 1.0)
    model.settings.temperature_field = openmc.TemperatureField(mesh, [294.0])

    with openmc.lib.TemporarySession(model, output=False, args=['-c']):
        field = openmc.lib.temperature_field
        assert field.size == 1
        assert len(field) == 1
        assert field[0] == 294.0
        field[0] = 294.0
        assert field[0] == 294.0
