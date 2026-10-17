from openmc.executor import _get_error_message
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
    # Any other failure
    ('Segmentation fault', 'OpenMC aborted unexpectedly.'),
])
def test_error_message_from_output(output, message):
    """The error message is recovered from each way OpenMC can fail."""
    assert _get_error_message(output) == message
