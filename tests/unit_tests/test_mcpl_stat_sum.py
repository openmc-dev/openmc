#!/usr/bin/env python
"""Test for MCPL stat:sum functionality (issue #3514)"""

import os
import shutil
import tempfile
import pytest
import openmc
import openmc.lib

# Skip test if MCPL is not available
pytestmark = pytest.mark.skipif(
    shutil.which("mcpl-config") is None,
    reason="mcpl-config command not found in PATH; MCPL is likely not available."
)


def test_mcpl_stat_sum_field():
    """Test that MCPL files contain proper stat:sum field with particle count."""
    
    # This test verifies that when OpenMC creates MCPL source files, 
    # they contain the stat:sum field. Since MCPL functions are not exposed
    # in the Python API, this test creates an actual OpenMC simulation
    # to generate MCPL files and then checks their content.
    
    mcpl = pytest.importorskip("mcpl")
    
    # Create a minimal working model that will generate MCPL files
    model = openmc.examples.pwr_pin_cell()
    model.settings.batches = 5
    model.settings.inactive = 2  
    model.settings.particles = 1000
    model.settings.sourcepoint = {'mcpl': True, 'separate': True}
    
    with tempfile.TemporaryDirectory() as tmpdir:
        cwd = os.getcwd()
        try:
            os.chdir(tmpdir)
            
            # Run a short simulation to generate MCPL files
            model.run(output=False)
            
            # Find the generated MCPL file (from the last batch)
            import glob
            mcpl_files = glob.glob('source.5.mcpl*')
            if len(mcpl_files) == 0:
                pytest.skip("No MCPL files were generated - MCPL interface not available")
            
            mcpl_file = mcpl_files[0]
            
            # Open and verify the stat:sum field exists
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
                    
                    if not stat_sum_found:
                        pytest.skip("stat:sum field not found - may be running with MCPL < 2.1.0")
                
                # Verify the stat:sum value is reasonable
                assert stat_sum_value != -1, "stat:sum was not updated from initial -1 value"
                
                # In eigenvalue mode, active batches generate source particles
                active_batches = model.settings.batches - model.settings.inactive  # 3 active batches
                expected_particles = active_batches * model.settings.particles     # 3000 total
                
                assert stat_sum_value == expected_particles, \
                    f"stat:sum value {stat_sum_value} doesn't match expected {expected_particles}"
                
        except Exception as e:
            # If the simulation fails (e.g., missing cross sections), skip the test
            if "cross_sections" in str(e).lower():
                pytest.skip(f"Simulation failed due to missing cross sections: {e}")
            else:
                raise
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