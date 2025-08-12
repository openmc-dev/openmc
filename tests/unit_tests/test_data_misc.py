#!/usr/bin/env python

from math import log
import os
from pathlib import Path

import numpy as np
import pytest
import openmc.data


def test_data_library(tmpdir):
    lib = openmc.data.DataLibrary.from_xml()
    for f in lib:
        assert sorted(f.keys()) == ['materials', 'path', 'type']

    f = lib.get_by_material('U235')
    assert f['type'] == 'neutron'
    assert 'U235' in f['materials']

    f = lib.get_by_material('c_H_in_H2O', data_type='thermal')
    assert f['type'] == 'thermal'
    assert 'c_H_in_H2O' in f['materials']

    lib.remove_by_material('Pu239')
    assert lib.get_by_material('Pu239') is None

    filename = str(tmpdir.join('test.xml'))
    lib.export_to_xml(filename)
    assert os.path.exists(filename)

    new_lib = openmc.data.DataLibrary()
    directory = os.path.dirname(os.environ['OPENMC_CROSS_SECTIONS'])
    new_lib.register_file(os.path.join(directory, 'H1.h5'))
    assert new_lib[-1]['type'] == 'neutron'
    new_lib.register_file(os.path.join(directory, 'c_Zr_in_ZrH.h5'))
    assert new_lib[-1]['type'] == 'thermal'


def test_depletion_chain_data_library(run_in_tmpdir):
    dep_lib = openmc.data.DataLibrary.from_xml()
    prev_len = len(dep_lib.libraries)
    chain_path = Path(__file__).parents[1] / "chain_simple.xml"
    dep_lib.register_file(chain_path)
    assert len(dep_lib.libraries) == prev_len + 1
    # Inspect
    dep_dict = dep_lib.libraries[-1]
    assert dep_dict['materials'] == []
    assert dep_dict['type'] == 'depletion_chain'
    assert dep_dict['path'] == str(chain_path)

    out_path = "cross_section_chain.xml"
    dep_lib.export_to_xml(out_path)

    dep_import = openmc.data.DataLibrary.from_xml(out_path)
    for lib in reversed(dep_import.libraries):
        if lib['type'] == 'depletion_chain':
            break
    else:
        raise IndexError("depletion_chain not found in exported DataLibrary")
    assert os.path.exists(lib['path'])


def test_linearize():
    """Test linearization of a continuous function."""
    x, y = openmc.data.linearize([-1., 1.], lambda x: 1 - x*x)
    f = openmc.data.Tabulated1D(x, y)
    assert f(-0.5) == pytest.approx(1 - 0.5*0.5, 0.001)
    assert f(0.32) == pytest.approx(1 - 0.32*0.32, 0.001)


def test_thin():
    """Test thinning of a tabulated function."""
    x = np.linspace(0., 2*np.pi, 1000)
    y = np.sin(x)
    x_thin, y_thin = openmc.data.thin(x, y)
    f = openmc.data.Tabulated1D(x_thin, y_thin)
    assert f(1.0) == pytest.approx(np.sin(1.0), 0.001)


def test_atomic_mass():
    assert openmc.data.atomic_mass('H1') == 1.007825031898
    assert openmc.data.atomic_mass('U235') == 235.043928117
    assert openmc.data.atomic_mass('Li6') == 6.01512288742
    assert openmc.data.atomic_mass('Pb220') == 220.025905
    with pytest.raises(KeyError):
        openmc.data.atomic_mass('U100')


def test_atomic_weight():
    assert openmc.data.atomic_weight('C') == 12.011115164865895
    assert openmc.data.atomic_weight('carbon') == 12.011115164865895
    with pytest.raises(ValueError):
        openmc.data.atomic_weight('Qt')


def test_water_density():
    dens = openmc.data.water_density
    # These test values are from IAPWS R7-97(2012).  They are actually specific
    # volumes so they need to be inverted.  They also need to be divided by 1000
    # to convert from [kg / m^3] to [g / cm^3].
    assert dens(300.0, 3.0) == pytest.approx(1e-3/0.100215168e-2, 1e-6)
    assert dens(300.0, 80.0) == pytest.approx(1e-3/0.971180894e-3, 1e-6)
    assert dens(500.0, 3.0) == pytest.approx(1e-3/0.120241800e-2, 1e-6)


def test_gnds_name():
    assert openmc.data.gnds_name(1, 1) == 'H1'
    assert openmc.data.gnds_name(40, 90) == ('Zr90')
    assert openmc.data.gnds_name(95, 242, 0) == ('Am242')
    assert openmc.data.gnds_name(95, 242, 1) == ('Am242_m1')
    assert openmc.data.gnds_name(95, 242, 10) == ('Am242_m10')


def test_isotopes():
    hydrogen_isotopes = [('H1', 0.99984426), ('H2', 0.00015574)]
    assert openmc.data.isotopes('H') == hydrogen_isotopes
    assert openmc.data.isotopes('hydrogen') == hydrogen_isotopes
    assert openmc.data.isotopes('Al') == [('Al27', 1.0)]
    assert openmc.data.isotopes('Aluminum') == [('Al27', 1.0)]
    assert openmc.data.isotopes('aluminium') == [('Al27', 1.0)]
    with pytest.raises(ValueError):
        openmc.data.isotopes('Чорнобиль')


def test_zam():
    assert openmc.data.zam('H1') == (1, 1, 0)
    assert openmc.data.zam('Zr90') == (40, 90, 0)
    assert openmc.data.zam('Am242') == (95, 242, 0)
    assert openmc.data.zam('Am242_m1') == (95, 242, 1)
    assert openmc.data.zam('Am242_m10') == (95, 242, 10)
    with pytest.raises(ValueError):
        openmc.data.zam('garbage')
    with pytest.raises(ValueError):
        openmc.data.zam('Am242-m1')

def test_half_life():
    with openmc.config.patch('chain_file', None):
        assert openmc.data.half_life('H2') is None
        assert openmc.data.half_life('U235') == pytest.approx(2.22102e16)
        assert openmc.data.half_life('Am242') == pytest.approx(57672.0)
        assert openmc.data.half_life('Am242_m1') == pytest.approx(4449622000.0)
        assert openmc.data.decay_constant('H2') == 0.0
        assert openmc.data.decay_constant('U235') == pytest.approx(log(2.0)/2.22102e16)
        assert openmc.data.decay_constant('Am242') == pytest.approx(log(2.0)/57672.0)
    assert openmc.data.decay_constant('Am242_m1') == pytest.approx(log(2.0)/4449622000.0)
