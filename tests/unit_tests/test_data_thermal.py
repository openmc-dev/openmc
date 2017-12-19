#!/usr/bin/env python

from collections import Callable
import os

import numpy as np
import pytest
import openmc.data


_ENDF_DATA = os.environ['OPENMC_ENDF_DATA']


@pytest.fixture(scope='module')
def h2o():
    """H in H2O thermal scattering data."""
    directory = os.path.dirname(os.environ['OPENMC_CROSS_SECTIONS'])
    filename = os.path.join(directory, 'c_H_in_H2O.h5')
    return openmc.data.ThermalScattering.from_hdf5(filename)


@pytest.fixture(scope='module')
def graphite():
    """Graphite thermal scattering data."""
    directory = os.path.dirname(os.environ['OPENMC_CROSS_SECTIONS'])
    filename = os.path.join(directory, 'c_Graphite.h5')
    return openmc.data.ThermalScattering.from_hdf5(filename)


@pytest.fixture(scope='module')
def h2o_njoy():
    path_h1 = os.path.join(_ENDF_DATA, 'neutrons', 'n-001_H_001.endf')
    path_h2o = os.path.join(_ENDF_DATA, 'thermal_scatt', 'tsl-HinH2O.endf')
    return openmc.data.ThermalScattering.from_njoy(
        path_h1, path_h2o, temperatures=[293.6, 500.0])


def test_h2o_attributes(h2o):
    assert h2o.name == 'c_H_in_H2O'
    assert h2o.nuclides == ['H1']
    assert h2o.secondary_mode == 'skewed'
    assert h2o.temperatures == ['294K']
    assert h2o.atomic_weight_ratio == pytest.approx(0.999167)


def test_h2o_xs(h2o):
    assert not h2o.elastic_xs
    for temperature, func in h2o.inelastic_xs.items():
        assert temperature.endswith('K')
        assert isinstance(func, Callable)


def test_graphite_attributes(graphite):
    assert graphite.name == 'c_Graphite'
    assert graphite.nuclides == ['C0', 'C12', 'C13']
    assert graphite.secondary_mode == 'skewed'
    assert graphite.temperatures == ['296K']
    assert graphite.atomic_weight_ratio == pytest.approx(11.898)


def test_graphite_xs(graphite):
    for temperature, func in graphite.elastic_xs.items():
        assert temperature.endswith('K')
        assert isinstance(func, openmc.data.CoherentElastic)
    for temperature, func in graphite.inelastic_xs.items():
        assert temperature.endswith('K')
        assert isinstance(func, Callable)
    elastic = graphite.elastic_xs['296K']
    assert elastic([1e-3, 1.0]) == pytest.approx([13.47464936, 0.62590156])


def test_export_to_hdf5(tmpdir, h2o_njoy, graphite):
    filename = str(tmpdir.join('water.h5'))
    h2o_njoy.export_to_hdf5(filename)
    assert os.path.exists(filename)

    filename = str(tmpdir.join('graphite.h5'))
    graphite.export_to_hdf5(filename)
    assert os.path.exists(filename)


def test_continuous_dist(h2o_njoy):
    for temperature, dist in h2o_njoy.inelastic_dist.items():
        assert temperature.endswith('K')
        assert isinstance(dist, openmc.data.CorrelatedAngleEnergy)


def test_get_thermal_name():
    f = openmc.data.get_thermal_name
    # Names which are recognized
    assert f('lwtr') == 'c_H_in_H2O'
    assert f('hh2o') == 'c_H_in_H2O'

    with pytest.warns(UserWarning, match='is not recognized'):
        # Names which can be guessed
        assert f('lw00') == 'c_H_in_H2O'
        assert f('graphite') == 'c_Graphite'
        assert f('D_in_D2O') == 'c_D_in_D2O'

        # Names that don't remotely match anything
        assert f('boogie_monster') == 'c_boogie_monster'
