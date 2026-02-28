from pytest import approx, raises

from openmc.data import mass_energy_absorption_coefficient
from openmc.data.function import Tabulated1D


def test_mass_energy_absorption_coefficient():
    # Spot checks on values from NIST tables
    mu_en = mass_energy_absorption_coefficient("air")
    assert isinstance(mu_en, Tabulated1D)
    assert mu_en.x[0] == approx(1e3)
    assert mu_en.y[0] == approx(3.599e3)
    assert mu_en.x[-1] == approx(2e7)
    assert mu_en.y[-1] == approx(1.311e-2)

    # Invalid material/data_source should raise an exception
    with raises(ValueError):
        mass_energy_absorption_coefficient("pasta")
    with raises(ValueError):
        mass_energy_absorption_coefficient("air", data_source="nist000")
