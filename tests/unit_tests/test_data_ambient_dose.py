from openmc.data import ambient_dose_coefficients
from pytest import approx, raises


def test_dose_coefficients():
    # Spot checks on values from ICRP tables

    energy, dose = ambient_dose_coefficients('photon', data_source='icrp74')
    assert energy[0] == approx(0.01e6)
    assert dose[0] == approx(0.061)
    assert energy[-1] == approx(10.0e6)
    assert dose[-1] == approx(25.6)

    energy, dose = ambient_dose_coefficients('neutron', data_source='icrp74')
    assert energy[0] == approx(1e-3)
    assert dose[0] == approx(6.60)
    assert energy[-1] == approx(20.0e6)
    assert dose[-1] == approx(600)

    # Invalid particle/data source should raise an exception
    with raises(ValueError):
        ambient_dose_coefficients('slime', 'icrp74')
    with raises(ValueError):
        ambient_dose_coefficients('neutron', data_source='icrp7000')
