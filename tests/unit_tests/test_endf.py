from openmc.data import endf
from pytest import approx


def test_float_endf():
    assert endf.float_endf('+3.2146') == approx(3.2146)
    assert endf.float_endf('.12345') == approx(0.12345)
    assert endf.float_endf('6.022+23') == approx(6.022e23)
    assert endf.float_endf('6.022-23') == approx(6.022e-23)
    assert endf.float_endf(' +1.01+ 2') == approx(101.0)
    assert endf.float_endf(' -1.01- 2') == approx(-0.0101)


def test_int_endf():
    assert endf.int_endf('    ') == 0
    assert endf.int_endf('+4032') == 4032
