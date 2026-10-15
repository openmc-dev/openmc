from pathlib import Path
import sys

import openmc
import pytest


@pytest.mark.parametrize('output, message', [
    # OpenMC's own error messages
    (' ERROR: Could not add nuclide', 'Could not add nuclide'),
    # libstdc++ (e.g. Linux)
    ("terminate called after throwing an instance of 'std::runtime_error'\n"
     "  what():  Could not add nuclide", 'Could not add nuclide'),
    # libc++ (e.g. macOS)
    ('libc++abi: terminating due to uncaught exception of type '
     'std::runtime_error: Could not add nuclide', 'Could not add nuclide'),
    # libc++ with an exception that has no message
    ('libc++abi: terminating due to uncaught exception of type int',
     'OpenMC aborted unexpectedly.'),
])
def test_error_message_from_output(run_in_tmpdir, output, message):
    """The error message is recovered from each way OpenMC can fail."""
    # Stand-in for the openmc executable that prints the output and fails
    fake_openmc = Path.cwd() / 'fake_openmc'
    with fake_openmc.open('w') as f:
        f.write(f'#!{sys.executable}\nimport sys\nprint({output!r})\nsys.exit(1)\n')
    fake_openmc.chmod(0o755)

    with pytest.raises(RuntimeError) as excinfo:
        openmc.plot_geometry(output=False, openmc_exec=str(fake_openmc))
    assert str(excinfo.value) == message
