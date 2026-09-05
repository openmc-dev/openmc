"""Tests for resuming an eigenvalue simulation from a statepoint."""

import re

import numpy as np
import pytest

import openmc


# One line is printed per generation, e.g.
#       9/1    1.18342    1.18055 +/- 0.00287
# capturing the batch, the generation and the running average k-effective.
_GENERATION_LINE = re.compile(
    r"^\s*(\d+)/(\d+)\s+[\d.]+\s+([\d.]+)\s+\+/-\s+[\d.]+\s*$", re.MULTILINE
)


@pytest.mark.parametrize("generations_per_batch", [1, 5])
def test_restart_keff(run_in_tmpdir, capsys, generations_per_batch):
    """A restarted run must resume with the correct accumulated k-effective.

    ``simulation::k_generation`` holds one entry per generation rather than one
    per batch, so the accumulation in ``restart_set_keff()`` has to be indexed
    by generation. Indexing it by batch sums the wrong entries -- and too few of
    them -- leaving the restarted run with a source-normalization k-effective
    roughly a factor of ``generations_per_batch`` too low, which over-produces
    fission sites and fills the shared fission bank.
    """
    model = openmc.examples.pwr_pin_cell()
    model.settings.particles = 500
    model.settings.inactive = 3
    model.settings.batches = 8
    model.settings.generations_per_batch = generations_per_batch

    statepoint = model.run()

    with openmc.StatePoint(statepoint) as sp:
        assert sp.generations_per_batch == generations_per_batch
        n_inactive = sp.n_inactive
        n_loaded = len(sp.k_generation)

    # Restart from the statepoint, asking for two more batches
    model.settings.batches = 10
    model.export_to_model_xml()
    capsys.readouterr()
    openmc.run(restart_file=statepoint)
    output = capsys.readouterr().out

    # Symptom of a low restart k-effective, asserted separately so that the
    # failure mode is obvious
    assert "fission bank is full" not in output

    with openmc.StatePoint("statepoint.10.h5") as sp:
        k_generation = np.asarray(sp.k_generation)

    average_k = [float(m[2]) for m in _GENERATION_LINE.findall(output)]
    assert len(average_k) == 2 * generations_per_batch

    # The running average reported for each generation of the restarted run must
    # be the mean over every active generation simulated so far, including the
    # ones loaded from the statepoint. The comparison is exact to the precision
    # of the printed value.
    first_active = generations_per_batch * n_inactive
    for i, reported in enumerate(average_k):
        expected = k_generation[first_active:n_loaded + i + 1].mean()
        assert reported == pytest.approx(expected, abs=1e-5)
