from pytest import approx, raises

from openmc.data import mass_energy_absorption_coefficient, mass_attenuation_coefficient
from openmc.data.function import Tabulated1D


def test_mass_attenuation_type():
    mu = mass_attenuation_coefficient(26)  # Fe
    assert isinstance(mu, Tabulated1D)


def test_mass_attenuation_spot_values():
    # Spot checks for Fe (Z=26) against NIST data: first/last tabulated points
    # and a mid-range value at 1 MeV
    mu = mass_attenuation_coefficient(26)
    assert mu(1e3) == approx(9085.0)
    assert mu(1e6) == approx(0.05995)
    assert mu(2e7) == approx(0.03224)


def test_mass_attenuation_caching():
    # Repeated calls with the same Z should return the identical object
    mu1 = mass_attenuation_coefficient(26)
    mu2 = mass_attenuation_coefficient('Fe')
    assert mu1 is mu2


def test_mass_attenuation_invalid_z():
    with raises(ValueError, match="Z=0"):
        mass_attenuation_coefficient(0)
    with raises(ValueError, match="Z=200"):
        mass_attenuation_coefficient(200)


def test_mass_energy_absorption_type():
    # Spot checks on values from NIST tables
    mu_en = mass_energy_absorption_coefficient("air")
    assert isinstance(mu_en, Tabulated1D)


def test_mass_energy_absorption_spot_values():
    mu_en = mass_energy_absorption_coefficient("air")
    assert mu_en(1e3) == approx(3.599e3)
    assert mu_en(10.e3) == approx(4.742)
    assert mu_en(2e7) == approx(1.311e-2)


def test_mass_energy_absorption_invalid():
    # Invalid material/data_source should raise an exception
    with raises(ValueError):
        mass_energy_absorption_coefficient("pasta")
    with raises(ValueError):
        mass_energy_absorption_coefficient("air", data_source="nist000")
