from pytest import approx, raises

from openmc.data import mu_en_coefficients


def test_mu_en_coefficients():
    # Spot checks on values from NIST tables
    energy, mu_en = mu_en_coefficients("air")
    assert energy[0] == approx(1e3)
    assert mu_en[0] == approx(3.599e3)
    assert energy[-1] == approx(2e7)
    assert mu_en[-1] == approx(1.311e-2)

    # Invalid particle/geometry should raise an exception
    with raises(ValueError):
        mu_en_coefficients("pasta")
    with raises(ValueError):
        mu_en_coefficients("air", data_source="nist000")
