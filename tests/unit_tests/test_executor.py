import os
import sys

import openmc
import pytest


@pytest.mark.parametrize('output', [
    # OpenMC's own error messages
    ' ERROR: Could not add nuclide',
    # libstdc++ (e.g. Linux)
    "terminate called after throwing an instance of 'std::runtime_error'\n"
    "  what():  Could not add nuclide",
    # libc++ (e.g. macOS)
    'libc++abi: terminating due to uncaught exception of type '
    'std::runtime_error: Could not add nuclide',
])
def test_error_message_from_output(run_in_tmpdir, output):
    """The error message is recovered from each way OpenMC can fail."""
    # Stand-in for the openmc executable that prints the output and fails
    fake_openmc = os.path.abspath('fake_openmc')
    with open(fake_openmc, 'w') as f:
        f.write(f'#!{sys.executable}\nimport sys\nprint({output!r})\nsys.exit(1)\n')
    os.chmod(fake_openmc, 0o755)

    with pytest.raises(RuntimeError, match='^Could not add nuclide$'):
        openmc.plot_geometry(output=False, openmc_exec=fake_openmc)
