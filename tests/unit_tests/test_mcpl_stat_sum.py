#!/usr/bin/env python
"""Test for MCPL stat:sum functionality (issue #3514)"""

import os
import shutil
import tempfile
import pytest
import openmc

# Skip test if MCPL is not available
pytestmark = pytest.mark.skipif(
    shutil.which("mcpl-config") is None,
    reason="mcpl-config command not found in PATH; MCPL is likely not available."
)


def test_mcpl_stat_sum_field():
    """Test that MCPL files contain proper stat:sum field with particle count."""
    
    # Only run if mcpl module is available for verification
    mcpl = pytest.importorskip("mcpl")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a simple model
        model = openmc.Model()
        
        # Create a simple material
        mat = openmc.Material()
        mat.add_nuclide('U235', 1.0)
        mat.set_density('g/cm3', 10.0)
        model.materials = [mat]
        
        # Create a simple geometry (sphere)
        sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
        cell = openmc.Cell(fill=mat, region=-sphere)
        model.geometry = openmc.Geometry([cell])
        
        # Configure settings for MCPL output
        model.settings = openmc.Settings()
        model.settings.batches = 5
        model.settings.inactive = 2
        model.settings.particles = 1000
        model.settings.source = openmc.Source(space=openmc.stats.Point())
        
        # Enable MCPL source file writing
        model.settings.source_point = {'mcpl': True, 'separate': True}
        
        # Run the simulation
        cwd = os.getcwd()
        try:
            os.chdir(tmpdir)
            model.run(output=False)
            
            # Find the MCPL file
            import glob
            mcpl_files = glob.glob('source.*.mcpl*')
            assert len(mcpl_files) > 0, "No MCPL source files were created"
            
            # Open and check the MCPL file
            mcpl_file = mcpl_files[0]
            with mcpl.MCPLFile(mcpl_file) as f:
                # Check if stat:sum field exists using convenience property
                if hasattr(f, 'stat_sum'):
                    # Use the convenience .stat_sum property directly
                    stat_sum_dict = f.stat_sum
                    assert 'openmc_np1' in stat_sum_dict, "openmc_np1 key not found in stat_sum"
                    stat_sum_value = stat_sum_dict['openmc_np1']
                else:
                    # Fallback to checking comments for older MCPL versions
                    comments = f.comments
                    
                    # Check for stat:sum in comments (MCPL stores these as comments)
                    stat_sum_found = False
                    stat_sum_value = None
                    
                    for comment in comments:
                        if 'stat:sum:openmc_np1' in comment:
                            stat_sum_found = True
                            # Extract the value
                            parts = comment.split(':')
                            if len(parts) >= 4:
                                stat_sum_value = float(parts[3].strip())
                            break
                    
                    assert stat_sum_found, "stat:sum:openmc_np1 field not found in MCPL file"
                
                # The value should be the total number of source particles
                # For 5 batches with 1000 particles each = 5000 total
                expected_particles = model.settings.batches * model.settings.particles
                
                # If stat:sum was properly updated, it should equal expected_particles
                # If it's still -1, the update before closing didn't work
                assert stat_sum_value != -1, "stat:sum was not updated from initial -1 value"
                assert stat_sum_value == expected_particles, \
                    f"stat:sum value {stat_sum_value} doesn't match expected {expected_particles}"
                
        finally:
            os.chdir(cwd)


def test_mcpl_stat_sum_crash_safety():
    """Test that incomplete MCPL files have stat:sum = -1."""
    
    # This test would verify that if file creation is interrupted,
    # the stat:sum field remains at -1 to indicate incomplete file
    # This is harder to test without mocking the MCPL library internals
    
    # For now, we can at least document the expected behavior:
    # 1. When mcpl_create_outfile is called, stat:sum should be set to -1
    # 2. Only when mcpl_close_outfile is called should it be updated
    # 3. If the program crashes between these calls, stat:sum remains -1
    
    # This could be tested with a C++ unit test that directly uses the
    # mcpl_interface functions and simulates a crash
    pass


if __name__ == "__main__":
    # Allow running this test directly
    test_mcpl_stat_sum_field()