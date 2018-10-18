#!/usr/bin/env python

from collections.abc import Mapping, Callable
import os

import numpy as np
import pandas as pd
import pytest
import openmc.data


_ENDF_DATA = os.environ['OPENMC_ENDF_DATA']


@pytest.fixture(scope='module')
def elements_endf():
    """Dictionary of element ENDF data indexed by atomic symbol."""
    elements = {'H': 1, 'O': 8, 'Al': 13, 'Cu': 29, 'Ag': 47, 'U': 92, 'Pu': 94}
    data = {}
    for symbol, Z in elements.items():
        p_file = 'photoat-{:03}_{}_000.endf'.format(Z, symbol)
        p_path = os.path.join(_ENDF_DATA, 'photoat', p_file)
        a_file = 'atom-{:03}_{}_000.endf'.format(Z, symbol)
        a_path = os.path.join(_ENDF_DATA, 'atomic_relax', a_file)
        data[symbol] = openmc.data.IncidentPhoton.from_endf(p_path, a_path)
    return data


@pytest.fixture()
def element(request, elements_endf):
    """Element ENDF data"""
    return elements_endf[request.param]


@pytest.mark.parametrize(
    'element, atomic_number', [
        ('Al', 13),
        ('Cu', 29),
        ('Pu', 94)
    ],
    indirect=['element']
)
def test_attributes(element, atomic_number):
    assert element.atomic_number == atomic_number


@pytest.mark.parametrize(
    'element, subshell, binding_energy, num_electrons', [
        ('H', 'K', 13.61, 1.0),
        ('O', 'L3', 14.15, 2.67),
        ('U', 'P2', 34.09, 2.0)
    ],
    indirect=['element']
)
def test_atomic_relaxation(element, subshell, binding_energy, num_electrons):
    atom_relax = element.atomic_relaxation
    assert isinstance(atom_relax, openmc.data.photon.AtomicRelaxation)
    assert subshell in atom_relax.subshells
    assert atom_relax.binding_energy[subshell] == binding_energy
    assert atom_relax.num_electrons[subshell] == num_electrons


@pytest.mark.parametrize('element', ['Al', 'Cu', 'Pu'], indirect=True)
def test_transitions(element):
    transitions = element.atomic_relaxation.transitions
    assert transitions
    assert isinstance(transitions, Mapping)
    for matrix in transitions.values():
        assert isinstance(matrix, pd.core.frame.DataFrame)
        assert len(matrix.columns) == 4
        assert sum(matrix['probability']) == pytest.approx(1.0)


@pytest.mark.parametrize('element', ['H', 'Al', 'Ag'], indirect=True)
def test_bremsstrahlung(element):
    brems = element.bremsstrahlung
    assert isinstance(brems, Mapping)
    assert np.all(np.diff(brems['electron_energy']) > 0.0)
    assert np.all(np.diff(brems['photon_energy']) > 0.0)
    assert brems['photon_energy'][0] == 0.0
    assert brems['photon_energy'][-1] == 1.0
    assert brems['dcs'].shape == (200, 30)


@pytest.mark.parametrize(
    'element, n_shell', [
        ('H', 1),
        ('O', 3),
        ('Al', 5)
    ],
    indirect=['element']
)
def test_compton_profiles(element, n_shell):
    profile = element.compton_profiles
    assert profile
    assert isinstance(profile, Mapping)
    assert all(isinstance(x, Callable) for x in profile['J'])
    assert all(len(x) == n_shell for x in profile.values())


@pytest.mark.parametrize(
    'element, reaction', [
        ('Cu', 541),
        ('Ag', 502),
        ('Pu', 504)
    ],
    indirect=['element']
)
def test_reactions(element, reaction):
    reactions = element.reactions
    assert all(isinstance(x, openmc.data.PhotonReaction) for x in reactions.values())
    assert reaction in reactions
    with pytest.raises(KeyError):
        reactions[18]


@pytest.mark.parametrize(
    'element, I', [
        ('H', 19.2),
        ('O', 95.0),
        ('U', 890.0)
    ],
    indirect=['element']
)
def test_stopping_powers(element, I):
    stopping_powers = element.stopping_powers
    assert isinstance(stopping_powers, Mapping)
    assert stopping_powers['I'] == I
    assert np.all(np.diff(stopping_powers['energy']) > 0.0)
    assert len(stopping_powers['s_collision']) == 200
    assert len(stopping_powers['s_radiative']) == 200


@pytest.mark.parametrize('element', ['Pu'], indirect=True)
def test_export_to_hdf5(tmpdir, element):
    filename = str(tmpdir.join('tmp.h5'))
    element.export_to_hdf5(filename)
    assert os.path.exists(filename)
