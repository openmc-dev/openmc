from pytest import approx, raises

from openmc.data import mu_en_coefficients


def test_mu_en_coefficients():
    # Spot checks on values from NIST tables
    energy, mu_en = mu_en_coefficients("air")
    assert energy[0] == approx(1e3)
    assert mu_en[0] == approx(3.599e3)
    assert energy[-1] == approx(2e7)
    assert mu_en[-1] == approx(1.311e-2)

    energy, mu_en = mu_en_coefficients("water")
    assert energy[0] == approx(1e3)
    assert mu_en[0] == approx(4.065e03)
    assert energy[-1] == approx(2e7)
    assert mu_en[-1] == approx(1.382e-2)

    energy, mu_en = mu_en_coefficients("water", data_source="nist126")
    assert energy[2] == approx(2e3)
    assert mu_en[2] == approx(6.152e02)
    assert energy[-2] == approx(1.5e7)
    assert mu_en[-2] == approx(1.441e-2)

    # Invalid particle/geometry should raise an exception
    with raises(ValueError) as excinfo:
        mu_en_coefficients("pasta")
    expected_materials = [
        "air",
        "water",
    ]
    expected_msg = (
        f"Unable to set 'material' to 'pasta' since it is not in {expected_materials}"
    )
    assert str(excinfo.value) == expected_msg
    with raises(ValueError) as excinfo:
        mu_en_coefficients("air", data_source="nist000")
    expected_msg = (
        f"Unable to set 'data_source' to 'nist000' since it is not in '{'nist126'}'"
    )
    assert str(excinfo.value) == expected_msg
