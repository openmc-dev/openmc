import numpy as np
import pytest

import openmc
import openmc.stats


def make_source(**kwargs):
    """Build a valid TokamakSource, overriding defaults via kwargs."""
    r_over_a = np.linspace(0.0, 1.0, 10)
    params = dict(
        major_radius=620.0,
        minor_radius=200.0,
        elongation=1.8,
        triangularity=0.45,
        shafranov_shift=10.0,
        r_over_a=r_over_a,
        emission_density=(1.0 - r_over_a**2),
        energy=openmc.stats.muir(e0=14.08e6, m_rat=5.0, kt=2.0e4),
    )
    params.update(kwargs)
    return openmc.TokamakSource(**params)


def test_tokamak_source_roundtrip():
    src = make_source(
        phi_start=0.1, phi_extent=np.pi, n_alpha=51, vertical_shift=5.0,
        strength=2.0, time=openmc.stats.Uniform(0.0, 1e-6))

    elem = src.to_xml_element()
    assert elem.get('type') == 'tokamak'

    new = openmc.SourceBase.from_xml_element(elem)
    assert isinstance(new, openmc.TokamakSource)
    assert new.major_radius == src.major_radius
    assert new.minor_radius == src.minor_radius
    assert new.elongation == src.elongation
    assert new.triangularity == src.triangularity
    assert new.shafranov_shift == src.shafranov_shift
    assert new.phi_start == src.phi_start
    assert new.phi_extent == src.phi_extent
    assert new.n_alpha == src.n_alpha
    assert new.vertical_shift == src.vertical_shift
    assert new.strength == src.strength
    np.testing.assert_allclose(new.r_over_a, src.r_over_a)
    np.testing.assert_allclose(new.emission_density, src.emission_density)
    assert len(new.energy) == 1
    assert isinstance(new.time, openmc.stats.Uniform)
    assert new.time.a == src.time.a
    assert new.time.b == src.time.b


def test_tokamak_source_default_time():
    src = make_source()
    assert src.time is None

    new = openmc.SourceBase.from_xml_element(src.to_xml_element())
    assert new.time is None

    with pytest.raises(TypeError):
        make_source(time=1.0)


def test_tokamak_source_multiple_energies():
    r_over_a = np.linspace(0.0, 1.0, 5)
    energies = [openmc.stats.muir(e0=14.08e6, m_rat=5.0, kt=kt)
                for kt in (1.0e4, 1.5e4, 2.0e4, 2.5e4, 3.0e4)]
    src = make_source(r_over_a=r_over_a,
                      emission_density=np.ones_like(r_over_a),
                      energy=energies)
    assert len(src.energy) == len(r_over_a)

    new = openmc.SourceBase.from_xml_element(src.to_xml_element())
    assert len(new.energy) == len(r_over_a)


@pytest.mark.parametrize("kwargs, match", [
    (dict(minor_radius=700.0), "smaller than major_radius"),
    (dict(shafranov_shift=150.0), "half the minor_radius"),
    (dict(emission_density=np.ones(5)), "same length as r_over_a"),
    (dict(energy=[openmc.stats.muir(14.08e6, 5.0, 2.0e4)] * 2),
     "Number of energy distributions"),
    (dict(r_over_a=np.linspace(0.1, 1.0, 10)), "must start at 0"),
    (dict(r_over_a=np.linspace(0.0, 0.9, 10)), "must end at 1"),
    (dict(emission_density=-np.linspace(0.0, 1.0, 10)), "cannot be negative"),
])
def test_tokamak_source_invalid(kwargs, match):
    with pytest.raises(ValueError, match=match):
        make_source(**kwargs)


@pytest.mark.parametrize("value", [-1.5, 1.5])
def test_tokamak_source_invalid_triangularity(value):
    with pytest.raises(ValueError):
        make_source(triangularity=value)


def test_tokamak_source_invalid_n_alpha():
    with pytest.raises(ValueError):
        make_source(n_alpha=2)
