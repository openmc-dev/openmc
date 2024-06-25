from openmc.data import dose_coefficients
from pytest import approx, raises


def test_dose_coefficients():
    # Spot checks on values from ICRP tables
    energy, dose = dose_coefficients('photon', 'AP', 'icrp116')
    assert energy[0] == approx(0.01e6)
    assert dose[0] == approx(0.0685)
    assert energy[-1] == approx(10e9)
    assert dose[-1] == approx(90.4)  # updated in corrigendum

    energy, dose = dose_coefficients('neutron', 'LLAT', 'icrp116')
    assert energy[0] == approx(1e-3)
    assert dose[0] == approx(1.04)
    assert energy[-1] == approx(10e9)
    assert dose[-1] == approx(1.23e3)

    energy, dose = dose_coefficients('neutron', 'LLAT', 'icrp119')
    assert energy[0] == approx(1e-3)
    assert dose[0] == approx(1.36)
    assert energy[-1] == approx(180e6)
    assert dose[-1] == approx(542.0)

    # the ISO column in icrp119 has NaN values, this test checks they are removed
    energy, dose = dose_coefficients('neutron', 'ISO', 'icrp119')
    assert energy[-1] == approx(20e6)
    assert dose[-1] == approx(343)
    assert len(energy) == len(dose)

    energy, dose = dose_coefficients('electron', 'ISO', 'icrp116')
    assert energy[0] == approx(0.01e6)
    assert dose[0] == approx(0.0188)
    assert energy[-1] == approx(10e9)
    assert dose[-1] == approx(699.0)

    # Invalid particle/geometry should raise an exception
    with raises(ValueError):
        dose_coefficients('slime', 'LAT')
    with raises(ValueError):
        dose_coefficients('neutron', 'ZZ')
    with raises(ValueError):
        dose_coefficients('neutron', 'ISO', 'foo')
