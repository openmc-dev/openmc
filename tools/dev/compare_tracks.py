#!/usr/bin/env python3
"""Report the first exact difference between two OpenMC track files."""

import argparse
import json
from pathlib import Path
import re

import h5py
import numpy as np
import openmc


# Logical order, independent of compound-member offsets and HDF5 byte order.
_FIELDS = (
    'r.x', 'r.y', 'r.z', 'u.x', 'u.y', 'u.z', 'E', 'time', 'wgt',
    'cell_id', 'cell_instance', 'material_id',
)
_FLOAT_FIELDS = _FIELDS[:9]
_SIGN = 1 << 63
_EXPONENT = 0x7ff0000000000000


def _field(array, name):
    for component in name.split('.'):
        array = array[component]
    return array


def _validate_dtype(dtype, fields):
    """Validate typed members, not padding, offsets, or byte order."""
    if dtype.names is None or set(dtype.names) != set(fields):
        raise ValueError(f'Expected state fields {tuple(fields)}')
    for name, expected in fields.items():
        member = dtype.fields[name][0]
        if isinstance(expected, dict):
            _validate_dtype(member, expected)
        elif member.subdtype or (member.kind, member.itemsize) != expected:
            raise ValueError(f'Invalid state dtype for {name}: {member}')


def _load_tracks(path):
    """Validate structural metadata before the Tracks reader slices records."""
    xyz = dict.fromkeys(('x', 'y', 'z'), ('f', 8))
    fields = {'r': xyz, 'u': xyz, 'E': ('f', 8), 'time': ('f', 8),
              'wgt': ('f', 8), 'cell_id': ('i', 4),
              'cell_instance': ('i', 4), 'material_id': ('i', 4)}
    identifiers = set()
    with h5py.File(path, 'r') as fh:
        version = np.asarray(fh.attrs['version'])
        if (version.ndim != 1 or version.size != 2
                or version.dtype.kind not in 'iu'
                or tuple(version) not in ((3, 0), (3, 1))):
            raise ValueError('Expected track format version 3.0 or 3.1')
        openmc.checkvalue.check_filetype_version(fh, 'track', 3)
        for name in fh:
            match = re.fullmatch(r'track_([0-9]+)_([0-9]+)_([0-9]+)', name)
            if match is None:
                raise ValueError(f'Invalid history name: {name}')
            identifier = tuple(map(int, match.groups()))
            if identifier in identifiers:
                raise ValueError(f'Duplicate history identifier: {identifier}')
            identifiers.add(identifier)
            if not isinstance(fh.get(name, getlink=True), h5py.HardLink):
                raise ValueError(f'{name}: linked track data are not supported')
            dset = fh[name]
            if not isinstance(dset, h5py.Dataset) or dset.ndim != 1:
                raise ValueError(f'{name}: expected a one-dimensional dataset')
            if dset.is_virtual or dset.external:
                raise ValueError(f'{name}: expected self-contained track data')
            _validate_dtype(dset.dtype, fields)
            particles = np.asarray(dset.attrs['particles'])
            offsets = np.asarray(dset.attrs['offsets'])
            count = np.asarray(dset.attrs['n_particles'])
            if (count.ndim != 0 or count.dtype.kind not in 'iu'
                    or particles.ndim != 1 or particles.dtype.kind not in 'iu'
                    or int(count) != particles.size or particles.size == 0):
                raise ValueError(f'{name}: inconsistent particle count')
            # Native tracks include at least their final state. Comparing
            # adjacent entries also avoids unsigned subtraction overflow.
            if (offsets.ndim != 1 or offsets.dtype.kind not in 'iu'
                    or offsets.size != particles.size + 1
                    or offsets[0] != 0 or offsets[-1] != len(dset)
                    or np.any(offsets[1:] <= offsets[:-1])):
                raise ValueError(f'{name}: invalid particle offsets')
            # Track 3.0 used legacy indices; 3.1 uses PDG numbers. The Tracks
            # reader normalizes the former to PDG. Validate the encoding first
            # so that this normalization cannot hide malformed 3.1 metadata.
            legacy = np.isin(particles, [0, 1, 2, 3])
            if version[1] == 0 and not np.all(legacy):
                raise ValueError(f'{name}: expected legacy particle indices')
            if version[1] == 1 and np.any(legacy):
                raise ValueError(f'{name}: expected PDG particle numbers')
    return openmc.Tracks(path)


def _bits(values):
    """View binary64 bits without floating arithmetic or NaN canonicalization."""
    return values.view(np.dtype('u8').newbyteorder(values.dtype.byteorder))


def _value(value, floating):
    if not floating:
        return {'value': int(value)}
    bits = int(_bits(np.asarray(value)))
    return {'value': repr(float(value)), 'bits': f'0x{bits:016x}'}


def _ulp_distance(left, right):
    """Count finite representable steps, with the two signed zeros coalesced."""
    if left & _EXPONENT == _EXPONENT or right & _EXPONENT == _EXPONENT:
        return None
    left_rank = -(left & (_SIGN - 1)) if left & _SIGN else left
    right_rank = -(right & (_SIGN - 1)) if right & _SIGN else right
    return abs(left_rank - right_rank)


def _first_state_difference(left, right):
    """Find the earliest state, breaking ties in documented field order."""
    count = min(len(left), len(right))
    first = None
    for order, name in enumerate(_FIELDS):
        a, b = _field(left[:count], name), _field(right[:count], name)
        floating = name in _FLOAT_FIELDS
        aa, bb = (_bits(a), _bits(b)) if floating else (a, b)
        different = np.flatnonzero(aa != bb)
        if different.size:
            index = int(different[0])
            if first is None or (index, order) < first[:2]:
                detail = {'kind': 'state', 'state_index': index, 'field': name,
                          'left': _value(a[index], floating),
                          'right': _value(b[index], floating)}
                if floating:
                    detail['ulp_distance'] = _ulp_distance(
                        int(aa[index]), int(bb[index]))
                first = (index, order, detail)
    return None if first is None else first[2]


def _summary(tracks):
    return {'histories': len(tracks),
            'particles': sum(len(t) for t in tracks),
            'states': sum(len(p.states) for t in tracks for p in t)}


def compare_tracks(left, right):
    """Compare recorded histories, returning a JSON-serializable report.

    Parameters
    ----------
    left, right : str or pathlib.Path
        Self-contained version-3.0 or 3.1 track files. Both must be closed by the
        producing simulations before comparison.

    Returns
    -------
    dict
        ``status`` is ``match``, ``mismatch``, or ``invalid_input``. On mismatch,
        ``difference`` identifies the first difference in history identifier,
        particle-track sequence, state index, then field order. Secondary
        indices are recorded sequence positions, not stable genealogy IDs.

    """
    inputs = []
    for side, path in (('left', left), ('right', right)):
        try:
            inputs.append(_load_tracks(path))
        except (OSError, ValueError, TypeError, KeyError, IndexError,
                OverflowError) as exc:
            return {'status': 'invalid_input', 'side': side, 'error': str(exc)}

    a_tracks, b_tracks = inputs
    report = {'status': 'match', 'left': _summary(a_tracks),
              'right': _summary(b_tracks)}
    if not a_tracks and not b_tracks:
        # Empty files cannot establish agreement of any recorded history.
        return {'status': 'invalid_input',
                'error': 'Neither input contains a recorded history'}
    a_by_id = {t.identifier: t for t in a_tracks}
    b_by_id = {t.identifier: t for t in b_tracks}

    for identifier in sorted(a_by_id.keys() | b_by_id.keys()):
        a, b = a_by_id.get(identifier), b_by_id.get(identifier)
        location = {'history': list(identifier)}
        difference = None
        if a is None or b is None:
            difference = {'kind': 'missing_history',
                          'missing_from': 'left' if a is None else 'right'}
        else:
            for index, (pa, pb) in enumerate(zip(a, b)):
                location = {'history': list(identifier),
                            'particle_index': index}
                if pa.particle != pb.particle:
                    difference = {'kind': 'particle_type',
                                  'left': int(pa.particle),
                                  'right': int(pb.particle)}
                else:
                    location['particle_type'] = int(pa.particle)
                    difference = _first_state_difference(pa.states, pb.states)
                    if difference is None and len(pa.states) != len(pb.states):
                        difference = {'kind': 'state_count',
                                      'state_index': min(len(pa.states),
                                                         len(pb.states)),
                                      'left': len(pa.states),
                                      'right': len(pb.states)}
                if difference is not None:
                    break
            if difference is None and len(a) != len(b):
                location = {'history': list(identifier),
                            'particle_index': min(len(a), len(b))}
                difference = {'kind': 'particle_count',
                              'left': len(a), 'right': len(b)}
        if difference is not None:
            report['status'] = 'mismatch'
            report['difference'] = {**location, **difference}
            return report
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('left', type=Path)
    parser.add_argument('right', type=Path)
    args = parser.parse_args(argv)
    report = compare_tracks(args.left, args.right)
    print(json.dumps(report, indent=2, allow_nan=False))
    return {'match': 0, 'mismatch': 1, 'invalid_input': 2}[report['status']]


if __name__ == '__main__':
    raise SystemExit(main())
