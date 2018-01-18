#!/usr/bin/env python

from collections import Mapping
import os

import numpy as np
import pytest
import openmc.data


def test_data_library(tmpdir):
    lib = openmc.data.DataLibrary.from_xml()
    for f in lib.libraries:
        assert sorted(f.keys()) == ['materials', 'path', 'type']

    f = lib.get_by_material('U235')
    assert f['type'] == 'neutron'
    assert 'U235' in f['materials']

    f = lib.get_by_material('c_H_in_H2O')
    assert f['type'] == 'thermal'
    assert 'c_H_in_H2O' in f['materials']

    filename = str(tmpdir.join('test.xml'))
    lib.export_to_xml(filename)
    assert os.path.exists(filename)

    new_lib = openmc.data.DataLibrary()
    directory = os.path.dirname(os.environ['OPENMC_CROSS_SECTIONS'])
    new_lib.register_file(os.path.join(directory, 'H1.h5'))
    assert new_lib.libraries[-1]['type'] == 'neutron'
    new_lib.register_file(os.path.join(directory, 'c_Zr_in_ZrH.h5'))
    assert new_lib.libraries[-1]['type'] == 'thermal'


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
