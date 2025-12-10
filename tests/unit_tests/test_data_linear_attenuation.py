import os

import numpy as np
import pytest

import openmc.data
import openmc.data.photon_attenuation as linear_attenuation
import openmc.data.photon_attenuation as photon_att
from openmc.data import IncidentPhoton
from openmc.data.function import Sum
from openmc.data.library import DataLibrary
from openmc.data.photon_attenuation import linear_attenuation_xs
from openmc.exceptions import DataError

PHOTON_MTS = (502, 504, 515, 517, 522)


@pytest.fixture(scope="module")
def xs_filename():
    xs = os.environ.get("OPENMC_CROSS_SECTIONS")
    if xs is None:
        pytest.skip("OPENMC_CROSS_SECTIONS not set.")
    return xs


@pytest.fixture(scope="module")
def elements_photon_xs(xs_filename):
    """Dictionary of IncidentPhoton data indexed by atomic symbol."""
    lib = DataLibrary.from_xml(xs_filename)

    elements = ["H", "O", "Al", "C", "Ag", "U", "Pb", "V"]
    data = {}
    for symbol in elements:
        entry = lib.get_by_material(symbol, data_type="photon")
        if entry is None:
            continue
        data[symbol] = IncidentPhoton.from_hdf5(entry["path"])
    return data


@pytest.mark.parametrize("symbol", ["C", "Pb"])
def test_linear_attenuation_xs_matches_sum(elements_photon_xs, symbol, monkeypatch):
    """linear_attenuation_xs should reproduce the sum of the relevant
    reaction channels from IncidentPhoton.reactions.
    """
    element = elements_photon_xs.get(symbol)
    if element is None:
        pytest.skip(f"No photon data for {symbol} in cross section library.")

    assert isinstance(element, openmc.data.IncidentPhoton)

    # Use preloaded IncidentPhoton instead of reading via DataLibrary in the helper
    monkeypatch.setattr(linear_attenuation, "_get_photon_data", lambda _: element)

    xs_sum = linear_attenuation_xs(symbol, temperature=293.6)

    # If the element has no relevant reactions, helper should return None
    has_relevant = any(mt in element.reactions for mt in PHOTON_MTS)
    if not has_relevant:
        assert xs_sum is None
        return

    assert isinstance(xs_sum, Sum)

    # Compare against explicit sum of reaction cross sections
    energy = np.logspace(2, 4, 50)
    expected = np.zeros_like(energy)
    for mt in PHOTON_MTS:
        if mt in element.reactions:
            expected += element.reactions[mt].xs(energy)

    actual = xs_sum(energy)
    assert np.allclose(actual, expected)


def test_linear_attenuation_xs_returns_none_when_no_photon_data(monkeypatch):
    """If _get_photon_data returns None, the helper should return None."""
    monkeypatch.setattr(linear_attenuation, "_get_photon_data", lambda _: None)

    xs_sum = linear_attenuation_xs("NonExistent", temperature=300.0)
    assert xs_sum is None


# ================================================================
# Tests for _get_photon_data (internal helper)
# ================================================================


def test_get_photon_data_valid(xs_filename):
    """_get_photon_data should load an IncidentPhoton object from the
    cross sections library and cache it.
    """
    lib = DataLibrary.from_xml(xs_filename)

    photon_nuclides = [mat for mat in lib if "photon" in mat["type"]]
    if not photon_nuclides:
        pytest.skip("No photon data entries available in cross section library.")

    nuclide = photon_nuclides[0]["materials"][0]

    # Clear internal cache
    photon_att._PHOTON_LIB = None
    photon_att._PHOTON_DATA = {}

    # Call target function
    data1 = photon_att._get_photon_data(nuclide)

    assert isinstance(data1, IncidentPhoton)

    # Cached instance should be reused on repeated calls
    data2 = photon_att._get_photon_data(nuclide)
    assert data1 is data2  # same object, cached


def test_get_photon_data_missing_nuclide():
    """_get_photon_data should return None when the nuclide has no photon data."""
    photon_att._PHOTON_LIB = None
    photon_att._PHOTON_DATA = {}

    # Pick a nuclide name guaranteed *not* to exist
    bad_name = "NonExistentNuclide_XXXX"

    data = photon_att._get_photon_data(bad_name)
    assert data is None


def test_get_photon_data_no_library(monkeypatch):
    """If DataLibrary.from_xml() fails, _get_photon_data should raise DataError."""
    # Force DataLibrary.from_xml to throw
    monkeypatch.setattr(
        photon_att.DataLibrary,
        "from_xml",
        lambda *_, **kw: (kw, (_ for _ in ()).throw(IOError("missing file")))[1],
    )

    # Clear caches
    photon_att._PHOTON_LIB = None
    photon_att._PHOTON_DATA = {}

    with pytest.raises(DataError):
        photon_att._get_photon_data("U235")


def test_linear_attenuation_reference_values(elements_photon_xs, monkeypatch):
    """Check linear_attenuation_xs for Pb and V at two reference energies."""
    pb_data = elements_photon_xs.get("Pb")
    v_data = elements_photon_xs.get("V")

    if pb_data is None or v_data is None:
        pytest.skip("Pb or V photon data not available in cross section library.")

    # Route _get_photon_data to our preloaded IncidentPhoton objects
    def _fake_get_photon_data(name: str):
        if name == "Pb":
            return pb_data
        if name == "V":
            return v_data
        return None

    monkeypatch.setattr(linear_attenuation, "_get_photon_data", _fake_get_photon_data)


    # Call the helper at room temperature
    xs_pb = linear_attenuation_xs("Pb", temperature=293.6)
    xs_v = linear_attenuation_xs("V", temperature=293.6)

    if xs_pb is None or xs_v is None:
        pytest.skip("No relevant photon reactions for Pb or V.")

    assert isinstance(xs_pb, Sum)
    assert isinstance(xs_v, Sum)

    # Test Lead
    pb_energies = np.array([1.0e5, 1.0e6])
    pb_vals = xs_pb(pb_energies)

    # data from https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z82.html
    expected_pb = np.array(
        [
            5.549e00,
            7.102e-02,
        ]
    )

    pb_mat = openmc.Material(temperature=293.6)
    pb_mat.add_element("Pb", 1.0)
    pb_mat.set_density("g/cm3", 11.34)

    expected_pb *= pb_mat.get_mass_density()/pb_mat.get_element_atom_densities()["Pb"]

    # Test Vanadium
    v_energies = np.array([1.0e5, 1.0e6])
    v_vals = xs_v(v_energies)

    # data from https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z23.html
    expected_v = np.array(
        [
            2.877e-01,
            5.794e-02,
        ]
    )

    v_mat = openmc.Material(temperature=293.6)
    v_mat.add_element("V", 1.0)
    v_mat.set_density("g/cm3", 11.34)

    expected_v *= pb_mat.get_mass_density()/v_mat.get_element_atom_densities()["V"]


    # Replace with tighter tolerances once real values are in
    assert np.allclose(pb_vals, expected_pb, rtol = 1e-4, atol=0)
    assert np.allclose(v_vals, expected_v, rtol = 1e-4, atol=0)
