"""
Test for the Random Ray Solver multiple point sources error detection fix.

This test verifies the fix for issue #3470 where multiple point sources 
placed in the same subdivided source region would silently overwrite each other.
The fix adds proper error detection to prevent this condition.
"""

import os
import tempfile
import pytest
import openmc
from openmc.examples import random_ray_three_region_cube


def test_multiple_point_sources_same_location_error(run_in_tmpdir):
    """Test that multiple point sources at the same location trigger an error 
    when mesh subdivision is enabled."""
    
    # Create a basic random ray model
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    
    # Enable mesh subdivision - this is required to trigger the error
    mesh = openmc.RegularMesh()
    mesh.dimension = (4, 4, 4)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (30.0, 30.0, 30.0)
    
    # Get the source universe to apply mesh subdivision
    source_universe = model.geometry.get_universes_by_name('source universe')[0]
    model.settings.random_ray['source_region_meshes'] = [
        (mesh, [source_universe])
    ]
    
    # Create multiple point sources at the exact same location
    # This should trigger the error when mesh subdivision is enabled
    point_location = openmc.stats.Point(xyz=(5.0, 5.0, 5.0))
    energy_dist = openmc.stats.Discrete([10.0], [1.0])
    
    source1 = openmc.IndependentSource(
        space=point_location, 
        energy=energy_dist, 
        strength=1.0
    )
    source2 = openmc.IndependentSource(
        space=point_location,  # Same location 
        energy=energy_dist, 
        strength=1.0
    )
    
    # Set the multiple point sources
    model.settings.source = [source1, source2]
    
    # Very short simulation to just trigger initialization
    model.settings.batches = 1
    model.settings.inactive = 0 
    model.settings.particles = 10
    
    # The error should be triggered during initialization/conversion of external sources
    with pytest.raises(RuntimeError, match=r"Multiple point sources detected.*same subdivided.*source region"):
        model.run()


def test_multiple_point_sources_different_locations_success(run_in_tmpdir):
    """Test that multiple point sources at different locations work correctly 
    with mesh subdivision enabled."""
    
    # Create a basic random ray model
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    
    # Enable mesh subdivision 
    mesh = openmc.RegularMesh()
    mesh.dimension = (4, 4, 4)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (30.0, 30.0, 30.0)
    
    # Get the source universe to apply mesh subdivision
    source_universe = model.geometry.get_universes_by_name('source universe')[0]
    model.settings.random_ray['source_region_meshes'] = [
        (mesh, [source_universe])
    ]
    
    # Create multiple point sources at different locations
    # This should work fine
    point_location1 = openmc.stats.Point(xyz=(2.0, 2.0, 2.0))
    point_location2 = openmc.stats.Point(xyz=(8.0, 8.0, 8.0))
    energy_dist = openmc.stats.Discrete([10.0], [1.0])
    
    source1 = openmc.IndependentSource(
        space=point_location1, 
        energy=energy_dist, 
        strength=1.0
    )
    source2 = openmc.IndependentSource(
        space=point_location2,  # Different location
        energy=energy_dist, 
        strength=1.0
    )
    
    # Set the multiple point sources
    model.settings.source = [source1, source2]
    
    # Very short simulation to just verify it initializes correctly
    model.settings.batches = 1
    model.settings.inactive = 0
    model.settings.particles = 10
    
    # This should not raise an error
    try:
        model.run()
    except RuntimeError as e:
        # If we get a RuntimeError, it should NOT be about multiple point sources
        assert "Multiple point sources detected" not in str(e)
        # Re-raise any other errors (could be due to test environment)
        raise


def test_multiple_point_sources_no_mesh_subdivision_success(run_in_tmpdir):
    """Test that multiple point sources at the same location work correctly 
    when mesh subdivision is disabled."""
    
    # Create a basic random ray model without mesh subdivision
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    
    # Explicitly ensure no mesh subdivision is configured
    # (default behavior, but being explicit for clarity)
    if 'source_region_meshes' in model.settings.random_ray:
        del model.settings.random_ray['source_region_meshes']
    
    # Create multiple point sources at the same location
    # This should work fine without mesh subdivision
    point_location = openmc.stats.Point(xyz=(5.0, 5.0, 5.0))
    energy_dist = openmc.stats.Discrete([10.0], [1.0])
    
    source1 = openmc.IndependentSource(
        space=point_location, 
        energy=energy_dist, 
        strength=1.0
    )
    source2 = openmc.IndependentSource(
        space=point_location,  # Same location, but no mesh subdivision
        energy=energy_dist, 
        strength=1.0
    )
    
    # Set the multiple point sources
    model.settings.source = [source1, source2]
    
    # Very short simulation to just verify it initializes correctly  
    model.settings.batches = 1
    model.settings.inactive = 0
    model.settings.particles = 10
    
    # This should not raise an error about multiple point sources
    try:
        model.run()
    except RuntimeError as e:
        # If we get a RuntimeError, it should NOT be about multiple point sources
        assert "Multiple point sources detected" not in str(e)
        # Re-raise any other errors (could be due to test environment)
        raise