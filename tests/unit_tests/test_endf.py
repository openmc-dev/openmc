from openmc.data import endf
from pytest import approx


def test_evaluation_from_material(endf_data):
    filename = f'{endf_data}/neutrons/n-001_H_001.endf'
    material = endf.get_evaluations(filename)[0]
    evaluation = endf.Evaluation(material)

    assert evaluation.material == material.MAT
    assert evaluation.gnds_name == 'H1'
    assert evaluation.section == material.section_text
    assert evaluation.reaction_list == material[1, 451]['section_list']


def test_float_endf():
    assert endf.float_endf('+3.2146') == approx(3.2146)
    assert endf.float_endf('.12345') == approx(0.12345)
    assert endf.float_endf('6.022+23') == approx(6.022e23)
    assert endf.float_endf('6.022-23') == approx(6.022e-23)
    assert endf.float_endf(' +1.01+ 2') == approx(101.0)
    assert endf.float_endf(' -1.01- 2') == approx(-0.0101)
    assert endf.float_endf('+ 2 . 3+ 1') == approx(23.0)
    assert endf.float_endf('-7 .8 -1') == approx(-0.78)
    assert endf.float_endf('3.14e0') == approx(3.14)
    assert endf.float_endf('3.14E0') == approx(3.14)
    assert endf.float_endf('3.14e-1') == approx(0.314)
    assert endf.float_endf('3.14d0') == approx(3.14)
    assert endf.float_endf('3.14D0') == approx(3.14)
    assert endf.float_endf('3.14d-1') == approx(0.314)
    assert endf.float_endf('1+2') == approx(100.0)
    assert endf.float_endf('-1+2') == approx(-100.0)
    assert endf.float_endf('1.+2') == approx(100.0)
    assert endf.float_endf('-1.+2') == approx(-100.0)
    assert endf.float_endf('        ') == 0.0
    assert endf.float_endf('9.876540000000000') == approx(9.87654)
    assert endf.float_endf('-2.225002+6') == approx(-2.225002e+6)


def test_int_endf():
    assert endf.int_endf('    ') == 0
    assert endf.int_endf('+4032') == 4032
