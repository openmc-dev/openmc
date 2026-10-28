from pathlib import Path

import openmc
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


@pytest.fixture
def tolerant_harness():
    return TolerantPyAPITestHarness('statepoint.1.h5', openmc.Model())


@pytest.mark.parametrize(('actual', 'expected', 'tolerance', 'accepted'), [
    ('inf', '1.0', 1e-6, False),
    ('inf', 'inf', 1e-6, True),
    ('nan', 'nan', 1e-6, False),
    ('1.0', '1.00001', 1e-6, False),
    ('1000000.0', '1000000.1', 1e-6, True),
    ('0.0', '-0.0', 1e-6, True),
])
def test_tolerant_harness_numeric_tokens(
        tolerant_harness, tmp_path, actual, expected, tolerance, accepted):
    """Non-finite values must not weaken the finite relative-tolerance gate."""
    actual_path = tmp_path / 'actual.dat'
    expected_path = tmp_path / 'expected.dat'
    actual_path.write_text(f'tally 1:\n{actual}\n')
    expected_path.write_text(f'tally 1:\n{expected}\n')

    assert tolerant_harness._are_files_equal(
        actual_path, expected_path, tolerance) is accepted


@pytest.mark.parametrize(('actual', 'expected', 'accepted'), [
    ('tally 1:\n1.0\n', 'tally 1:\n1.0\n', True),
    ('tally 1:\n1.0\n', 'tally 1:\n1.0\n2.0\n', False),
    ('tally 1:\n1.0 2.0\n', 'tally 1:\n1.0\n', False),
    ('tally 1:\n1.0\n', 'tally 1:\nnot-a-number\n', False),
])
def test_tolerant_harness_file_structure(
        tolerant_harness, tmp_path, actual, expected, accepted):
    """The numeric comparison must preserve text, line and token checks."""
    actual_path = tmp_path / 'actual.dat'
    expected_path = tmp_path / 'expected.dat'
    actual_path.write_text(actual)
    expected_path.write_text(expected)

    assert tolerant_harness._are_files_equal(
        actual_path, expected_path, 1e-6) is accepted


def test_tolerant_harness_rejects_and_preserves_results(
        tolerant_harness, run_in_tmpdir, capsys):
    """Rejected results must raise and retain the actual output for review."""
    actual_text = 'tally 1:\ninf\n'
    expected_text = 'tally 1:\n1.0\n'
    actual_path = Path('results_test.dat')
    expected_path = Path('results_true.dat')
    actual_path.write_text(actual_text)
    expected_path.write_text(expected_text)

    with pytest.raises(AssertionError, match='Results do not agree'):
        tolerant_harness._compare_results()

    assert not actual_path.exists()
    assert Path('results_error.dat').read_text() == actual_text
    assert expected_path.read_text() == expected_text
    assert 'Result differences:' in capsys.readouterr().out


def test_tolerant_harness_accepts_matching_infinity(
        tolerant_harness, run_in_tmpdir):
    """An infinity only matches an infinity with the same sign."""
    text = 'tally 1:\ninf\n'
    Path('results_test.dat').write_text(text)
    Path('results_true.dat').write_text(text)

    tolerant_harness._compare_results()

    assert Path('results_test.dat').read_text() == text
    assert not Path('results_error.dat').exists()
