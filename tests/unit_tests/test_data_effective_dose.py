from openmc.data import dose_coefficients
from pytest import approx, raises


def test_dose_coefficients():
    # Spot checks on values from ICRP tables
    energy, dose = dose_coefficients('photon', 'AP')
    assert energy[0] == approx(0.01e6)
    assert dose[0] == approx(0.0685)
    assert energy[-1] == approx(10e9)
    assert dose[-1] == approx(90.4)  # updated in corrigendum

    energy, dose = dose_coefficients('neutron', 'LLAT')
    assert energy[0] == approx(1e-3)
    assert dose[0] == approx(1.04)
    assert energy[-1] == approx(10e9)
    assert dose[-1] == approx(1.23e3)

    energy, dose = dose_coefficients('electron', 'ISO')
    assert energy[0] == approx(0.01e6)
    assert dose[0] == approx(0.0188)
    assert energy[-1] == approx(10e9)
    assert dose[-1] == approx(699.0)

    energy, dose = dose_coefficients('photon', data_source='icrp74')
    assert energy[0] == approx(0.01e6)
    assert dose[0] == approx(7.43*0.00653)
    assert energy[-1] == approx(10.0e6)
    assert dose[-1] == approx(24.0*0.990)

    energy, dose = dose_coefficients('neutron', 'LLAT', data_source='icrp74')
    assert energy[0] == approx(1e-3)
    assert dose[0] == approx(1.68)
    assert energy[-1] == approx(20.0e6)
    assert dose[-1] == approx(338.0)
    
    energy, dose = dose_coefficients('neutron', 'LLAT', data_source='icrp74')
    assert energy[0] == approx(1e-3)
    assert dose[0] == approx(1.68)
    assert energy[-1] == approx(20.0e6)
    assert dose[-1] == approx(338.0)
    
    energy, dose = dose_coefficients('neutron', 
                                     data_source='icrp74',
                                     dose_type = 'ambient')
    assert energy[0] == approx(1e-3)
    assert dose[0] == approx(6.60)
    assert energy[-1] == approx(20.0e6)
    assert dose[-1] == approx(600)
    
    energy, dose = dose_coefficients('photon', 
                                     data_source='icrp74',
                                     dose_type = 'ambient')
    assert energy[0] == approx(0.01e6)
    assert dose[0] == approx(0.061)
    assert energy[-1] == approx(10e6)
    assert dose[-1] == approx(25.6)

    # Invalid particle/geometry should raise an exception
    with raises(ValueError):
        dose_coefficients('slime', 'LAT')
    with raises(ValueError):
        dose_coefficients('neutron', 'ZZ')
    with raises(ValueError):
        dose_coefficients('neutron', data_source='icrp7000')
    with raises(ValueError):
        dose_coefficients('neutron', 
                          data_source='icrp116',
                          dose_type='ambient')
    with raises(ValueError) as excinfo:
        dose_coefficients("photons", data_source="icrp116")
    expected_particles = [
        "electron",
        "helium",
        "mu+",
        "mu-",
        "neutron",
        "photon",
        "photon kerma",
        "pi+",
        "pi-",
        "positron",
        "proton",
    ]
    expected_msg = (
        "'photons' has no effective dose data in data source icrp116."
        f" Available particles for icrp116 are: {expected_particles}"
    )
    assert str(excinfo.value) == expected_msg
