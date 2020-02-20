#!/usr/bin/env python

from collections.abc import Mapping, Callable
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import openmc.data


@pytest.fixture(scope='module')
def elements_endf():
    """Dictionary of element ENDF data indexed by atomic symbol."""
    endf_data = os.environ['OPENMC_ENDF_DATA']
    elements = {'H': 1, 'O': 8, 'Al': 13, 'Cu': 29, 'Ag': 47, 'U': 92, 'Pu': 94}
    data = {}
    for symbol, Z in elements.items():
        p_file = 'photoat-{:03}_{}_000.endf'.format(Z, symbol)
        p_path = os.path.join(endf_data, 'photoat', p_file)
        a_file = 'atom-{:03}_{}_000.endf'.format(Z, symbol)
        a_path = os.path.join(endf_data, 'atomic_relax', a_file)
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


@pytest.mark.parametrize(
    'element, I, i_shell, ionization_energy, num_electrons', [
        ('H', 19.2, 0, 13.6, 1),
        ('O', 95.0, 2, 13.62, 4),
        ('U', 890.0, 25, 6.033, -3)
    ],
    indirect=['element']
)
def test_bremsstrahlung(element, I, i_shell, ionization_energy, num_electrons):
    brems = element.bremsstrahlung
    assert isinstance(brems, Mapping)
    assert brems['I'] == I
    assert brems['num_electrons'][i_shell] == num_electrons
    assert brems['ionization_energy'][i_shell] == ionization_energy
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


@pytest.mark.parametrize('element', ['Pu'], indirect=True)
def test_export_to_hdf5(tmpdir, element):
    filename = str(tmpdir.join('tmp.h5'))
    element.export_to_hdf5(filename)
    assert os.path.exists(filename)
    # Read in data from hdf5
    element2 = openmc.data.IncidentPhoton.from_hdf5(filename)
    # Check for some cross section and datasets of element and element2
    energy = np.logspace(np.log10(1.0), np.log10(1.0e10), num=100)
    for mt in (502, 504, 515, 517, 522, 541, 570):
        xs = element[mt].xs(energy)
        xs2 = element2[mt].xs(energy)
        assert np.allclose(xs, xs2)
    assert element[502].scattering_factor == element2[502].scattering_factor
    assert element.atomic_relaxation.transitions['O3'].equals(
           element2.atomic_relaxation.transitions['O3'])
    assert (element.compton_profiles['binding_energy'] ==
           element2.compton_profiles['binding_energy']).all()
    assert (element.bremsstrahlung['electron_energy'] ==
           element2.bremsstrahlung['electron_energy']).all()
    # Export to hdf5 again
    element2.export_to_hdf5(filename, 'w')

def test_photodat_only(run_in_tmpdir):
    endf_dir = Path(os.environ['OPENMC_ENDF_DATA'])
    photoatomic_file = endf_dir / 'photoat' / 'photoat-001_H_000.endf'
    data = openmc.data.IncidentPhoton.from_endf(photoatomic_file)
    data.export_to_hdf5('tmp.h5', 'w')