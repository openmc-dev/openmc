from collections.abc import Mapping, Callable
from copy import deepcopy
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import openmc.data


@pytest.fixture(scope='module')
def elements_endf(endf_data):
    """Dictionary of element ENDF data indexed by atomic symbol."""
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


def test_photodat_only(run_in_tmpdir, endf_data):
    endf_dir = Path(endf_data)
    photoatomic_file = endf_dir / 'photoat' / 'photoat-001_H_000.endf'
    data = openmc.data.IncidentPhoton.from_endf(photoatomic_file)
    data.export_to_hdf5('tmp.h5', 'w')


def test_from_endf_material(endf_data):
    endf_dir = Path(endf_data)
    photoatomic_file = endf_dir / 'photoat' / 'photoat-001_H_000.endf'
    relaxation_file = endf_dir / 'atomic_relax' / 'atom-001_H_000.endf'
    photoatomic = openmc.data.endf.get_evaluations(photoatomic_file)[0]
    relaxation = openmc.data.endf.get_evaluations(relaxation_file)[0]

    data = openmc.data.IncidentPhoton.from_endf(photoatomic, relaxation)

    assert data.atomic_number == 1
    assert 502 in data.reactions
    assert data.atomic_relaxation.binding_energy['K'] == pytest.approx(13.61)


def test_atomic_relaxation_from_endf_material(endf_data):
    filename = Path(endf_data) / 'atomic_relax' / 'atom-001_H_000.endf'
    material = openmc.data.endf.get_evaluations(filename)[0]

    data = openmc.data.AtomicRelaxation.from_endf(material)

    assert data.binding_energy['K'] == pytest.approx(13.61)
    assert data.num_electrons['K'] == pytest.approx(1.0)


@pytest.fixture(scope='module')
def photon_evaluations(endf_data):
    endf_dir = Path(endf_data)
    paths = (
        endf_dir / 'photoat' / 'photoat-001_H_000.endf',
        endf_dir / 'atomic_relax' / 'atom-001_H_000.endf',
    )
    return tuple(openmc.data.endf.Evaluation(path) for path in paths)


@pytest.fixture
def photon_with_metadata(photon_evaluations):
    photoatomic, relaxation = deepcopy(photon_evaluations)
    photoatomic.info['library'] = ('Photoatomic evaluation', 8, 1)
    relaxation.info['library'] = ('Relaxation evaluation', 7, 3)
    return openmc.data.IncidentPhoton.from_endf(photoatomic, relaxation)


@pytest.mark.parametrize(
    'input_type', ['str', 'path', 'evaluation', 'material'])
def test_source_metadata_from_endf(endf_data, photon_evaluations, input_type):
    """Extract each component's own ENDF library, version and release."""
    inputs = [
        Path(endf_data) / 'photoat' / 'photoat-001_H_000.endf',
        Path(endf_data) / 'atomic_relax' / 'atom-001_H_000.endf',
    ]
    if input_type == 'str':
        inputs = [str(path) for path in inputs]
    elif input_type == 'evaluation':
        inputs = deepcopy(photon_evaluations)
    elif input_type == 'material':
        inputs = [openmc.data.endf.get_evaluations(path)[0] for path in inputs]

    data = openmc.data.IncidentPhoton.from_endf(*inputs)
    standalone = openmc.data.AtomicRelaxation.from_endf(inputs[1])
    for component, evaluation in zip(
            (data, data.atomic_relaxation), photon_evaluations):
        library, version, release = evaluation.info['library']
        assert component.source_metadata == {
            'library': library, 'version': version, 'release': release}
    assert standalone.source_metadata == data.atomic_relaxation.source_metadata


def test_source_metadata_components_are_independent(photon_evaluations):
    """Keep different evaluations and their mutable source records separate."""
    photoatomic, relaxation = deepcopy(photon_evaluations)
    photoatomic.info['library'] = ('Photoatomic evaluation', 8, 1)
    relaxation.info['library'] = ('Relaxation evaluation', 7, 3)
    data = openmc.data.IncidentPhoton.from_endf(photoatomic, relaxation)

    assert data.source_metadata == {
        'library': 'Photoatomic evaluation', 'version': 8, 'release': 1}
    assert data.atomic_relaxation.source_metadata == {
        'library': 'Relaxation evaluation', 'version': 7, 'release': 3}
    photoatomic.info['library'] = ('Changed input', 99, 99)
    data.source_metadata['version'] = 9
    assert data.source_metadata['library'] == 'Photoatomic evaluation'
    assert data.atomic_relaxation.source_metadata['version'] == 7
    assert relaxation.info['library'] == ('Relaxation evaluation', 7, 3)


def test_source_metadata_hdf5_roundtrip(tmp_path, photon_with_metadata):
    """Round-trip independent records through paths and HDF5 groups."""
    path = tmp_path / 'photon.h5'
    data = photon_with_metadata
    data.export_to_hdf5(path)

    for filename in (str(path), path):
        restored = openmc.data.IncidentPhoton.from_hdf5(filename)
        assert restored.source_metadata == data.source_metadata
        assert (restored.atomic_relaxation.source_metadata ==
                data.atomic_relaxation.source_metadata)
        for component in (restored, restored.atomic_relaxation):
            assert isinstance(component.source_metadata['library'], str)
            assert type(component.source_metadata['version']) is int
            assert type(component.source_metadata['release']) is int

    with h5py.File(path, 'r') as h5file:
        group = h5file[data.name]
        restored = openmc.data.IncidentPhoton.from_hdf5(group)
        standalone = openmc.data.AtomicRelaxation.from_hdf5(group['subshells'])
        assert restored.source_metadata == data.source_metadata
        assert (standalone.source_metadata ==
                data.atomic_relaxation.source_metadata)
        for location, metadata in (
                (group, data.source_metadata),
                (group['subshells'], data.atomic_relaxation.source_metadata)):
            library = location.attrs['source_library']
            if isinstance(library, bytes):
                library = library.decode('utf-8')
            assert library == metadata['library']
            assert location.attrs['source_version'] == metadata['version']
            assert location.attrs['source_release'] == metadata['release']
        assert h5file.id.valid


def test_source_metadata_legacy_hdf5(tmp_path, photon_with_metadata):
    """Read older files without inventing provenance for either source."""
    path = tmp_path / 'legacy.h5'
    data = photon_with_metadata
    data.export_to_hdf5(path)
    original = openmc.data.IncidentPhoton.from_hdf5(path)
    with h5py.File(path, 'r+') as h5file:
        for group in (h5file[data.name], h5file[data.name]['subshells']):
            for key in ('source_library', 'source_version', 'source_release'):
                del group.attrs[key]

    restored = openmc.data.IncidentPhoton.from_hdf5(path)
    assert restored.source_metadata == {}
    assert restored.atomic_relaxation.source_metadata == {}
    assert (restored.atomic_relaxation.binding_energy ==
            data.atomic_relaxation.binding_energy)
    np.testing.assert_array_equal(restored[502].xs.y, original[502].xs.y)


def test_source_metadata_preserves_hdf5_payload(
        tmp_path, photon_with_metadata):
    """Keep numerical data, attributes and element registration unchanged."""
    enriched = tmp_path / 'with_metadata.h5'
    plain = tmp_path / 'without_metadata.h5'
    data = photon_with_metadata
    data.export_to_hdf5(enriched)
    without_metadata = deepcopy(data)
    without_metadata.source_metadata = {}
    without_metadata.atomic_relaxation.source_metadata = {}
    without_metadata.export_to_hdf5(plain)

    metadata_keys = {'source_library', 'source_version', 'source_release'}
    with h5py.File(enriched, 'r') as actual, h5py.File(plain, 'r') as expected:
        actual_names, expected_names = [], []
        actual.visit(actual_names.append)
        expected.visit(expected_names.append)
        assert actual_names == expected_names
        assert list(actual) == [data.name]
        for name in ['', *actual_names]:
            left, right = actual[name or '/'], expected[name or '/']
            assert set(left.attrs) - metadata_keys == set(right.attrs)
            for key in right.attrs:
                np.testing.assert_array_equal(
                    left.attrs[key], right.attrs[key])
            if isinstance(right, h5py.Dataset):
                assert left.dtype == right.dtype
                np.testing.assert_array_equal(left[()], right[()])

    library = openmc.data.DataLibrary()
    library.register_file(enriched)
    assert len(library) == 1
    assert library[0]['type'] == 'photon'
    assert library[0]['materials'] == [data.name]


def test_source_metadata_without_relaxation(tmp_path, photon_evaluations):
    """Do not assign photoatomic provenance to absent relaxation data."""
    photoatomic, _ = photon_evaluations
    data = openmc.data.IncidentPhoton.from_endf(photoatomic)
    assert data.atomic_relaxation is None
    library, version, release = photoatomic.info['library']
    assert data.source_metadata == {
        'library': library, 'version': version, 'release': release}
    path = tmp_path / 'photoatomic.h5'
    data.export_to_hdf5(path)
    with h5py.File(path, 'r') as h5file:
        subshells = h5file[data.name]['subshells']
        assert not any(key.startswith('source_') for key in subshells.attrs)
    restored = openmc.data.IncidentPhoton.from_hdf5(path)
    assert restored.source_metadata == data.source_metadata
    assert restored.atomic_relaxation.source_metadata == {}


def test_source_metadata_defaults_are_independent():
    """Give new objects separate empty provenance dictionaries."""
    first = openmc.data.IncidentPhoton(1)
    second = openmc.data.IncidentPhoton(1)
    first_relaxation = openmc.data.AtomicRelaxation({}, {}, {})
    second_relaxation = openmc.data.AtomicRelaxation({}, {}, {})
    for component in (first, second, first_relaxation, second_relaxation):
        assert component.source_metadata == {}
    first.source_metadata['library'] = 'Photoatomic evaluation'
    first_relaxation.source_metadata['library'] = 'Relaxation evaluation'
    assert second.source_metadata == {}
    assert second_relaxation.source_metadata == {}


@pytest.mark.parametrize('as_bytes', [False, True])
def test_source_metadata_unicode(tmp_path, photon_with_metadata, as_bytes):
    """Decode UTF-8 provenance without changing the numeric version fields."""
    data = photon_with_metadata
    data.source_metadata['library'] = 'Évaluation photonique'
    path = tmp_path / 'unicode.h5'
    data.export_to_hdf5(path)
    if as_bytes:
        with h5py.File(path, 'r+') as h5file:
            attrs = h5file[data.name].attrs
            del attrs['source_library']
            attrs['source_library'] = np.bytes_(
                data.source_metadata['library'].encode('utf-8'))
    restored = openmc.data.IncidentPhoton.from_hdf5(path)
    assert restored.source_metadata == data.source_metadata


def test_source_metadata_partial_records(tmp_path, photon_with_metadata):
    """Preserve known fields without filling in absent component metadata."""
    data = photon_with_metadata
    data.source_metadata = {'library': 'Partial photoatomic record'}
    data.atomic_relaxation.source_metadata = {'version': 7}
    path = tmp_path / 'partial.h5'
    data.export_to_hdf5(path)
    restored = openmc.data.IncidentPhoton.from_hdf5(path)
    assert restored.source_metadata == data.source_metadata
    assert (restored.atomic_relaxation.source_metadata ==
            data.atomic_relaxation.source_metadata)


def test_source_metadata_zero_versions(tmp_path, photon_with_metadata):
    """Retain zero-valued NumPy integers as ordinary Python metadata."""
    data = photon_with_metadata
    data.source_metadata = {
        'library': 'Zero version',
        'version': np.int32(0),
        'release': np.int64(0),
    }
    path = tmp_path / 'zero.h5'
    data.export_to_hdf5(path)
    restored = openmc.data.IncidentPhoton.from_hdf5(path)
    assert restored.source_metadata == {
        'library': 'Zero version', 'version': 0, 'release': 0}
    assert type(restored.source_metadata['version']) is int
    assert type(restored.source_metadata['release']) is int


@pytest.mark.parametrize('component, metadata, error', [
    ('photoatomic', {'version': '8'}, TypeError),
    ('relaxation', {'library': 7}, TypeError),
    ('photoatomic', {'unrecognized': 'source'}, ValueError),
    ('photoatomic', {'library': 'invalid\x00library'}, ValueError),
    ('relaxation', {'library': '\ud800'}, UnicodeEncodeError),
    ('photoatomic', {'version': 1 << 100}, ValueError),
    ('relaxation', {'release': -(1 << 100)}, ValueError),
])
def test_invalid_source_metadata_preserves_output(
        tmp_path, photon_with_metadata, component, metadata, error):
    """Reject invalid metadata before truncating an existing output file."""
    data = photon_with_metadata
    target = data if component == 'photoatomic' else data.atomic_relaxation
    target.source_metadata = metadata
    path = tmp_path / 'existing.h5'
    original = b'Existing output must survive invalid metadata.'
    path.write_bytes(original)

    with pytest.raises(error):
        data.export_to_hdf5(path, mode='w')

    assert path.read_bytes() == original
