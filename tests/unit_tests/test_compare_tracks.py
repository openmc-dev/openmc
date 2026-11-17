"""Exact track diagnostics using self-contained, synthetic HDF5 records."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest


_SCRIPT = Path(__file__).parents[2] / 'tools/dev/compare_tracks.py'
_SPEC = importlib.util.spec_from_file_location('compare_tracks', _SCRIPT)
_TOOL = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_TOOL)

_FLOAT_FIELDS = (
    'r.x', 'r.y', 'r.z', 'u.x', 'u.y', 'u.z', 'E', 'time', 'wgt',
)
_INT_FIELDS = ('cell_id', 'cell_instance', 'material_id')


def _dtype(byteorder='<', padded=False, reversed_fields=False):
    xyz = [(axis, byteorder + 'f8') for axis in 'xyz']
    fields = [('r', xyz), ('u', xyz)]
    fields.extend((name, byteorder + 'f8') for name in ('E', 'time', 'wgt'))
    fields.extend((name, byteorder + 'i4') for name in _INT_FIELDS)
    if reversed_fields:
        fields.reverse()
    return np.dtype(fields, align=padded)


def _field(states, field):
    for name in field.split('.'):
        states = states[name]
    return states


def _states(size=3, **dtype_kwargs):
    states = np.zeros(size, dtype=_dtype(**dtype_kwargs))
    for index, field in enumerate(_FLOAT_FIELDS):
        _field(states, field)[:] = np.arange(size) + index + 0.5
    for index, field in enumerate(_INT_FIELDS):
        states[field] = np.arange(size) + index + 1
    return states


def _write(path, histories=None):
    if histories is None:
        histories = [('track_1_1_2', [(_states(), 2112)])]
    with h5py.File(path, 'w') as fh:
        fh.attrs['filetype'] = np.bytes_('track')
        fh.attrs['version'] = np.array([3, 1], dtype='i4')
        for name, particles in histories:
            states = np.concatenate(
                [state for state, _ in particles], dtype=particles[0][0].dtype)
            dset = fh.create_dataset(name, data=states)
            dset.attrs['n_particles'] = len(particles)
            dset.attrs['particles'] = [pdg for _, pdg in particles]
            dset.attrs['offsets'] = np.concatenate([
                [0], np.cumsum([len(state) for state, _ in particles]),
            ])
    return path


@pytest.fixture
def track_pair(tmp_path):
    return _write(tmp_path / 'left.h5'), _write(tmp_path / 'right.h5')


def _compare_states(tmp_path, left, right):
    a = _write(tmp_path / 'a.h5', [('track_1_1_2', [(left, 2112)])])
    b = _write(tmp_path / 'b.h5', [('track_1_1_2', [(right, 2112)])])
    return _TOOL.compare_tracks(a, b)


def test_match(track_pair):
    report = _TOOL.compare_tracks(*track_pair)
    assert report == {
        'status': 'match',
        'left': {'histories': 1, 'particles': 1, 'states': 3},
        'right': {'histories': 1, 'particles': 1, 'states': 3},
    }


@pytest.mark.parametrize('legacy,pdg', [(0, 2112), (1, 22), (2, 11), (3, -11)])
def test_legacy_and_current_particle_codes_match(tmp_path, legacy, pdg):
    a = _write(tmp_path / 'a.h5', [('track_1_1_2', [(_states(), legacy)])])
    b = _write(tmp_path / 'b.h5', [('track_1_1_2', [(_states(), pdg)])])
    with h5py.File(a, 'r+') as fh:
        fh.attrs['version'] = [3, 0]
    assert _TOOL.compare_tracks(a, b)['status'] == 'match'
    assert _TOOL.compare_tracks(b, a)['status'] == 'match'


def test_legacy_particle_type_difference(tmp_path):
    paths = []
    for legacy in (0, 1):
        path = _write(tmp_path / f'{legacy}.h5', [
            ('track_1_1_2', [(_states(), legacy)]),
        ])
        with h5py.File(path, 'r+') as fh:
            fh.attrs['version'] = [3, 0]
        paths.append(path)
    report = _TOOL.compare_tracks(*paths)
    assert report['status'] == 'mismatch'
    assert report['difference']['kind'] == 'particle_type'
    assert report['difference']['left'] == 2112
    assert report['difference']['right'] == 22


def test_legacy_file_rejects_mixed_particle_conventions(track_pair):
    _write(track_pair[1], [
        ('track_1_1_2', [(_states(), 0), (_states(), 22)]),
    ])
    with h5py.File(track_pair[1], 'r+') as fh:
        fh.attrs['version'] = [3, 0]
    report = _TOOL.compare_tracks(*track_pair)
    assert report['status'] == 'invalid_input'
    assert report['side'] == 'right'


@pytest.mark.parametrize('field', _FLOAT_FIELDS + _INT_FIELDS)
def test_each_state_field(tmp_path, field):
    left, right = _states(), _states()
    _field(right, field)[1] += 1
    report = _compare_states(tmp_path, left, right)
    assert report['status'] == 'mismatch'
    difference = report['difference']
    assert difference['history'] == [1, 1, 2]
    assert difference['particle_index'] == 0
    assert difference['particle_type'] == 2112
    assert difference['kind'] == 'state'
    assert difference['state_index'] == 1
    assert difference['field'] == field
    if field in _FLOAT_FIELDS:
        assert difference['left']['bits'].startswith('0x')
        assert difference['right']['bits'].startswith('0x')
        assert difference['ulp_distance'] > 0
    else:
        assert difference['left']['value'] + 1 == difference['right']['value']


def test_earliest_state_before_field(tmp_path):
    left, right = _states(), _states()
    right['r']['x'][2] += 1
    right['material_id'][0] += 1
    difference = _compare_states(tmp_path, left, right)['difference']
    assert difference['state_index'] == 0
    assert difference['field'] == 'material_id'


def test_same_state_field_order(tmp_path):
    left, right = _states(), _states()
    right['u']['y'][1] += 1
    right['u']['z'][1] += 1
    right['E'][1] += 1
    difference = _compare_states(tmp_path, left, right)['difference']
    assert (difference['state_index'], difference['field']) == (1, 'u.y')


@pytest.mark.parametrize('left_bits,right_bits,ulp', [
    (0x0000000000000000, 0x8000000000000000, 0),
    (0x0000000000000000, 0x0000000000000001, 1),
    (0x8000000000000001, 0x0000000000000001, 2),
    (0x0010000000000000, 0x000fffffffffffff, 1),
    (0x3ff0000000000000, 0x3ff0000000000001, 1),
    (0xbff0000000000000, 0xbff0000000000001, 1),
    (0xbff0000000000000, 0x3ff0000000000000,
     2 * 0x3ff0000000000000),
    (0x7fefffffffffffff, 0xffefffffffffffff,
     2 * 0x7fefffffffffffff),
    (0x7ff0000000000000, 0xfff0000000000000, None),
    (0x7fefffffffffffff, 0x7ff0000000000000, None),
    (0x7ff8000000000001, 0x7ff8000000000002, None),
    (0x7ff0000000000001, 0x7ff8000000000001, None),
    (0x7ff8000000000001, 0xfff8000000000001, None),
])
def test_exact_float_bits(tmp_path, left_bits, right_bits, ulp):
    left, right = _states(), _states()
    left['E'].view('<u8')[1] = left_bits
    right['E'].view('<u8')[1] = right_bits
    report = _compare_states(tmp_path, left, right)
    assert report['status'] == 'mismatch'
    difference = report['difference']
    assert difference['field'] == 'E'
    assert difference['left']['bits'] == f'0x{left_bits:016x}'
    assert difference['right']['bits'] == f'0x{right_bits:016x}'
    assert difference['ulp_distance'] == ulp
    # Nonfinite values must remain representable in standards-compliant JSON.
    assert json.loads(json.dumps(report, allow_nan=False)) == report


@pytest.mark.parametrize('bits', [
    0x8000000000000000, 0x7ff0000000000000,
    0xfff0000000000000, 0x7ff8000000000001, 0x7ff0000000000001,
])
def test_identical_nonfinite_and_signed_zero(tmp_path, bits):
    left, right = _states(), _states()
    left['E'].view('<u8')[1] = bits
    right['E'].view('<u8')[1] = bits
    assert _compare_states(tmp_path, left, right)['status'] == 'match'


def test_storage_byteorder_padding_and_member_order(tmp_path):
    left = _states()
    right = _states(byteorder='>', padded=True, reversed_fields=True)
    assert left.dtype != right.dtype
    assert left.dtype.itemsize != right.dtype.itemsize
    report = _compare_states(tmp_path, left, right)
    assert report['status'] == 'match'
    with h5py.File(tmp_path / 'b.h5', 'r') as fh:
        dtype = fh['track_1_1_2'].dtype
        assert dtype.names != left.dtype.names
        assert dtype['E'].byteorder == '>'
        assert dtype.itemsize != left.dtype.itemsize


def test_big_endian_nan_payload(tmp_path):
    left = _states()
    right = _states(byteorder='>')
    left['E'].view('<u8')[1] = 0x7ff8000000000001
    right['E'].view('>u8')[1] = 0x7ff8000000000001
    assert _compare_states(tmp_path, left, right)['status'] == 'match'
    right['E'].view('>u8')[1] = 0x7ff8000000000002
    difference = _compare_states(tmp_path, left, right)['difference']
    assert difference['right']['bits'] == '0x7ff8000000000002'


def test_numeric_history_order(tmp_path):
    left, right = _states(), _states()
    right['E'][0] += 1
    histories = [
        ('track_1_1_10', [(left, 2112)]),
        ('track_1_1_2', [(left, 2112)]),
    ]
    a = _write(tmp_path / 'a.h5', histories)
    b = _write(tmp_path / 'b.h5', [
        (name, [(right, 2112)]) for name, _ in reversed(histories)
    ])
    assert _TOOL.compare_tracks(a, b)['difference']['history'] == [1, 1, 2]


def test_particle_ordinal_before_state(tmp_path):
    left, right = _states(), _states()
    right['E'][2] += 1
    a = _write(tmp_path / 'a.h5', [
        ('track_1_1_2', [(left, 22), (left, 11)]),
    ])
    b = _write(tmp_path / 'b.h5', [
        ('track_1_1_2', [(right, 22), (right, 11)]),
    ])
    difference = _TOOL.compare_tracks(a, b)['difference']
    assert difference['particle_index'] == 0
    assert difference['particle_type'] == 22
    assert difference['state_index'] == 2


@pytest.mark.parametrize('missing_from', ['left', 'right'])
def test_missing_history(tmp_path, missing_from):
    a = _write(tmp_path / 'a.h5')
    b = _write(tmp_path / 'b.h5', [
        ('track_1_1_2', [(_states(), 2112)]),
        ('track_1_1_10', [(_states(), 2112)]),
    ])
    paths = (a, b) if missing_from == 'left' else (b, a)
    report = _TOOL.compare_tracks(*paths)
    assert report['status'] == 'mismatch'
    assert report['difference'] == {
        'history': [1, 1, 10], 'kind': 'missing_history',
        'missing_from': missing_from,
    }


@pytest.mark.parametrize('kind', ['particle_type', 'state_count',
                                 'particle_count'])
def test_structural_difference(tmp_path, kind):
    left = [(_states(), 2112), (_states(), 22)]
    right = [(_states(), 2112), (_states(), 22)]
    if kind == 'particle_type':
        right[1] = (_states(), -11)
    elif kind == 'state_count':
        right[1] = (_states(2), 22)
    else:
        right.pop()
    a = _write(tmp_path / 'a.h5', [('track_1_1_2', left)])
    b = _write(tmp_path / 'b.h5', [('track_1_1_2', right)])
    report = _TOOL.compare_tracks(a, b)
    assert report['status'] == 'mismatch'
    difference = report['difference']
    assert difference['kind'] == kind
    assert difference['particle_index'] == 1
    if kind == 'state_count':
        assert difference['state_index'] == 2
        assert (difference['left'], difference['right']) == (3, 2)
    elif kind == 'particle_count':
        assert (difference['left'], difference['right']) == (2, 1)
    else:
        assert (difference['left'], difference['right']) == (22, -11)
        assert 'particle_type' not in difference


def test_state_difference_before_count_difference(tmp_path):
    left, right = _states(), _states(2)
    right['time'][0] += 1
    difference = _compare_states(tmp_path, left, right)['difference']
    assert difference['kind'] == 'state'
    assert difference['state_index'] == 0
    assert difference['field'] == 'time'


@pytest.mark.parametrize('attribute,value', [
    ('offsets', [-1, 3]),
    ('offsets', [1, 3]),
    ('offsets', [0, 2]),
    ('offsets', [0, 4]),
    ('offsets', [0]),
    ('offsets', [0, 2, 3]),
    ('offsets', [0.0, 3.0]),
    ('offsets', [[0, 3]]),
    ('offsets', np.array([0, 2**63], dtype='u8')),
    ('particles', [0]),
    ('particles', [1]),
    ('particles', [2]),
    ('particles', [3]),
    ('particles', [2112.0]),
    ('particles', [[2112]]),
    ('particles', []),
    ('n_particles', -1),
    ('n_particles', 0),
    ('n_particles', 2),
    ('n_particles', [1]),
    ('n_particles', 1.0),
])
def test_invalid_particle_metadata(track_pair, attribute, value):
    with h5py.File(track_pair[1], 'r+') as fh:
        fh['track_1_1_2'].attrs[attribute] = value
    report = _TOOL.compare_tracks(*track_pair)
    assert report['status'] == 'invalid_input'
    assert report['side'] == 'right'


def test_nonmonotone_unsigned_offsets(track_pair):
    with h5py.File(track_pair[1], 'r+') as fh:
        attrs = fh['track_1_1_2'].attrs
        attrs['n_particles'] = 3
        attrs['particles'] = [2112, 22, 11]
        attrs['offsets'] = np.array([0, 2, 1, 3], dtype='u8')
    assert _TOOL.compare_tracks(*track_pair)['status'] == 'invalid_input'


@pytest.mark.parametrize('attribute', ['n_particles', 'particles', 'offsets'])
def test_missing_particle_metadata(track_pair, attribute):
    with h5py.File(track_pair[0], 'r+') as fh:
        del fh['track_1_1_2'].attrs[attribute]
    report = _TOOL.compare_tracks(*track_pair)
    assert report['status'] == 'invalid_input'
    assert report['side'] == 'left'


@pytest.mark.parametrize('attribute,value', [
    ('version', [2, 0]),
    ('version', [4, 0]),
    ('version', [3, 2]),
    ('version', [3, -1]),
    ('version', [3]),
    ('version', [3, 1, 0]),
    ('version', [3.0, 1.0]),
    ('version', []),
    ('version', [[3, 1]]),
    ('version', 3),
    ('filetype', np.bytes_('statepoint')),
])
def test_invalid_file_metadata(track_pair, attribute, value):
    with h5py.File(track_pair[1], 'r+') as fh:
        fh.attrs[attribute] = value
    assert _TOOL.compare_tracks(*track_pair)['status'] == 'invalid_input'


@pytest.mark.parametrize('attribute', ['version', 'filetype'])
def test_missing_file_metadata(track_pair, attribute):
    with h5py.File(track_pair[1], 'r+') as fh:
        del fh.attrs[attribute]
    assert _TOOL.compare_tracks(*track_pair)['status'] == 'invalid_input'


@pytest.mark.parametrize('malformation', [
    'group', 'scalar', 'matrix', 'wrong_float_width', 'wrong_integer_sign',
    'missing_field', 'extra_field', 'coordinate_array', 'bad_name',
    'duplicate_identifier', 'soft_link', 'external_link',
])
def test_invalid_dataset_schema(track_pair, malformation):
    with h5py.File(track_pair[1], 'r+') as fh:
        original = fh['track_1_1_2']
        states = original[()]
        attrs = dict(original.attrs)
        if malformation == 'bad_name':
            fh.move('track_1_1_2', 'track_one_1_2')
        elif malformation == 'duplicate_identifier':
            fh['track_01_1_2'] = original
        elif malformation == 'soft_link':
            fh['track_1_1_3'] = h5py.SoftLink('/track_1_1_2')
        elif malformation == 'external_link':
            fh['track_1_1_3'] = h5py.ExternalLink(
                str(track_pair[0]), '/track_1_1_2')
        else:
            del fh['track_1_1_2']
            if malformation == 'group':
                fh.create_group('track_1_1_2')
            elif malformation == 'scalar':
                fh.create_dataset('track_1_1_2', data=states[0])
            elif malformation == 'matrix':
                fh.create_dataset('track_1_1_2', data=states.reshape(1, 3))
            else:
                fields = states.dtype.descr
                if malformation == 'wrong_float_width':
                    fields = [(n, '<f4' if n == 'E' else t) for n, t in fields]
                elif malformation == 'wrong_integer_sign':
                    fields = [(n, '<u4' if n == 'cell_id' else t)
                              for n, t in fields]
                elif malformation == 'missing_field':
                    fields = [(n, t) for n, t in fields if n != 'time']
                elif malformation == 'extra_field':
                    fields.append(('extra', '<f8'))
                elif malformation == 'coordinate_array':
                    fields = [(n, ('<f8', (3,)) if n == 'r' else t)
                              for n, t in fields]
                dset = fh.create_dataset(
                    'track_1_1_2', data=np.zeros(3, dtype=fields))
                dset.attrs.update(attrs)
    assert _TOOL.compare_tracks(*track_pair)['status'] == 'invalid_input'


@pytest.mark.parametrize('empty_side', ['left', 'right', 'both'])
def test_empty_file_is_not_agreement(track_pair, empty_side):
    for index, side in enumerate(('left', 'right')):
        if empty_side in (side, 'both'):
            _write(track_pair[index], [])
    report = _TOOL.compare_tracks(*track_pair)
    if empty_side == 'both':
        assert report['status'] == 'invalid_input'
    else:
        assert report['status'] == 'mismatch'
        assert report['difference']['kind'] == 'missing_history'
        assert report['difference']['missing_from'] == empty_side


def test_late_malformed_history_takes_precedence(tmp_path):
    left, right = _states(), _states()
    right['E'][0] += 1
    a = _write(tmp_path / 'a.h5')
    b = _write(tmp_path / 'b.h5', [
        ('track_1_1_2', [(right, 2112)]),
        ('track_1_1_10', [(left, 2112)]),
    ])
    with h5py.File(b, 'r+') as fh:
        fh['track_1_1_10'].attrs['offsets'] = [-1, 3]
    report = _TOOL.compare_tracks(a, b)
    assert report['status'] == 'invalid_input'
    assert report['side'] == 'right'


def test_empty_particle_states_are_not_agreement(tmp_path):
    report = _compare_states(tmp_path, _states(0), _states(0))
    assert report['status'] == 'invalid_input'


@pytest.mark.parametrize('storage', ['external_link', 'external', 'virtual'])
def test_external_storage_rejected_before_read(track_pair, tmp_path, storage):
    with h5py.File(track_pair[1], 'r+') as fh:
        attrs = dict(fh['track_1_1_2'].attrs)
        del fh['track_1_1_2']
        absent = str(tmp_path / 'not-present.h5')
        if storage == 'external_link':
            fh['track_1_1_2'] = h5py.ExternalLink(absent, '/records')
        elif storage == 'external':
            dset = fh.create_dataset(
                'track_1_1_2', shape=(3,), dtype=_dtype(),
                external=[(absent, 0, h5py.h5f.UNLIMITED)],
            )
            dset.attrs.update(attrs)
        else:
            layout = h5py.VirtualLayout(shape=(3,), dtype=_dtype())
            layout[:] = h5py.VirtualSource(absent, 'records', shape=(3,))
            dset = fh.create_virtual_dataset('track_1_1_2', layout)
            dset.attrs.update(attrs)
    report = _TOOL.compare_tracks(*track_pair)
    assert report['status'] == 'invalid_input'
    expected = 'linked track data' if storage == 'external_link' else (
        'self-contained track data')
    assert expected in report['error']


def test_comparison_does_not_modify_inputs(track_pair):
    before = [path.read_bytes() for path in track_pair]
    assert _TOOL.compare_tracks(*track_pair)['status'] == 'match'
    assert [path.read_bytes() for path in track_pair] == before


def test_missing_file(track_pair, tmp_path):
    report = _TOOL.compare_tracks(track_pair[0], tmp_path / 'missing.h5')
    assert report['status'] == 'invalid_input'
    assert report['side'] == 'right'


@pytest.mark.parametrize('status,code', [
    ('match', 0), ('mismatch', 1), ('invalid_input', 2),
])
def test_cli_exit_status(track_pair, tmp_path, status, code):
    left, right = track_pair
    if status == 'mismatch':
        states = _states()
        states['E'][0] += 1
        _write(right, [('track_1_1_2', [(states, 2112)])])
    elif status == 'invalid_input':
        right = tmp_path / 'missing.h5'
    result = subprocess.run(
        [sys.executable, str(_SCRIPT), str(left), str(right)],
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == code, result.stderr
    assert json.loads(result.stdout)['status'] == status
    assert 'Traceback' not in result.stderr
