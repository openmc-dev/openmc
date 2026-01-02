"""Test for MCPL stat:sum functionality"""

from pathlib import Path
import shutil

import pytest
import openmc


@pytest.mark.skipif(shutil.which("mcpl-config") is None, reason="MCPL is not available.")
def test_mcpl_stat_sum_field(run_in_tmpdir):
    """Test that MCPL files contain proper stat:sum field with particle count.

    This test verifies that when OpenMC creates MCPL source files, they contain
    the stat:sum field. Since MCPL functions are not exposed in the Python API,
    this test creates an actual OpenMC simulation to generate MCPL files and
    then checks their content.
    """

    mcpl = pytest.importorskip("mcpl")

    # Create a minimal working model that will generate MCPL files
    model = openmc.examples.pwr_pin_cell()
    model.settings.batches = 5
    model.settings.inactive = 2
    model.settings.particles = 1000
    model.settings.sourcepoint = {'mcpl': True, 'separate': True}

    # Run a short simulation to generate MCPL files
    model.run(output=False)

    # Find the generated MCPL file (from the last batch)
    mcpl_file = Path('source.5.mcpl')
    assert mcpl_file.exists(), "No MCPL files were generated"

    # Open and verify the stat:sum field exists
    with mcpl.MCPLFile(mcpl_file) as f:
        # Check if stat:sum field exists using convenience property
        if hasattr(f, 'stat_sum'):
            # Use the convenience .stat_sum property directly
            stat_sum_dict = f.stat_sum
            assert 'openmc_np1' in stat_sum_dict, "openmc_np1 key not found in stat_sum"
            stat_sum_value = int(stat_sum_dict['openmc_np1'])
        else:
            # Fallback to checking comments for older MCPL versions
            comments = f.comments

            # Check for stat:sum in comments (MCPL stores these as comments)
            stat_sum_value = None

            for comment in comments:
                if 'stat:sum:openmc_np1' in comment:
                    # Extract the value
                    parts = comment.split(':')
                    if len(parts) >= 4:
                        stat_sum_value = int(parts[3].strip())
                    break
            else:
                pytest.skip("stat:sum field not found - may be running with MCPL < 2.1.0")

    # Verify the stat:sum value is reasonable
    assert stat_sum_value != -1, "stat:sum was not updated from initial -1 value"

    # In eigenvalue mode, active batches generate source particles
    active_batches = model.settings.batches - model.settings.inactive  # 3 active batches
    expected_particles = active_batches * model.settings.particles     # 3000 total

    assert stat_sum_value == expected_particles, \
        f"stat:sum value {stat_sum_value} doesn't match expected {expected_particles}"
