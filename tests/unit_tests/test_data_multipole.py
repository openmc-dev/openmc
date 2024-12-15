import os
import pathlib

import numpy as np
import pytest
import openmc.data


@pytest.fixture(scope="module")
def u235():
    directory = pathlib.Path(os.environ["OPENMC_CROSS_SECTIONS"]).parent
    u235 = directory / "wmp" / "092235.h5"
    return openmc.data.WindowedMultipole.from_hdf5(u235)


@pytest.fixture(scope="module")
def b10():
    directory = pathlib.Path(os.environ["OPENMC_CROSS_SECTIONS"]).parent
    b10 = directory / "wmp" / "005010.h5"
    return openmc.data.WindowedMultipole.from_hdf5(b10)


def test_evaluate(u235):
    """Test the cross section evaluation of a library."""
    energies = [1e-3, 1.0, 10.0, 50.0]
    scattering, absorption, fission = u235(energies, 0.0)
    assert (scattering[1], absorption[1], fission[1]) == pytest.approx(
        (13.09, 77.56, 67.36), rel=1e-3
    )
    scattering, absorption, fission = u235(energies, 300.0)
    assert (scattering[2], absorption[2], fission[2]) == pytest.approx(
        (11.24, 21.26, 15.50), rel=1e-3
    )


def test_evaluate_none_poles(b10):
    """Test a library with no poles, i.e., purely polynomials."""
    energies = [1e-3, 1.0, 10.0, 1e3, 1e5]
    scattering, absorption, fission = b10(energies, 0.0)
    assert (scattering[0], absorption[0], fission[0]) == pytest.approx(
        (2.201, 19330.0, 0.0), rel=1e-3
    )
    scattering, absorption, fission = b10(energies, 300.0)
    assert (scattering[-1], absorption[-1], fission[-1]) == pytest.approx(
        (2.878, 1.982, 0.0), rel=1e-3
    )


def test_export_to_hdf5(tmpdir, u235):
    filename = str(tmpdir.join("092235.h5"))
    u235.export_to_hdf5(filename)
    assert os.path.exists(filename)


def test_from_endf():
    pytest.importorskip("vectfit")
    endf_data = os.environ["OPENMC_ENDF_DATA"]
    endf_file = os.path.join(endf_data, "neutrons", "n-001_H_001.endf")
    return openmc.data.WindowedMultipole.from_endf(
        endf_file, log=True, wmp_options={"n_win": 400, "n_cf": 3}
    )


def test_from_endf_search():
    pytest.importorskip("vectfit")
    endf_data = os.environ["OPENMC_ENDF_DATA"]
    endf_file = os.path.join(endf_data, "neutrons", "n-095_Am_244.endf")
    return openmc.data.WindowedMultipole.from_endf(
        endf_file, log=True, wmp_options={"search": True, "rtol": 1e-2}
    )
