#!/usr/bin/env python3
"""Run a small simulation and print digests of its particle tracks and tallies.

The digests cover the raw bits of every recorded particle state, the k
effective estimates, and the tally results. Two platforms that produce the
same digests have produced bit-identical particle histories, which is the
property the portable math functions are meant to guarantee. When the digests
differ, tools such as tools/dev/compare_tracks.py (openmc-dev/openmc#4120) can
locate the first differing state.

Usage: track-digest.py OUTPUT.json [--particles N] [--batches N]
"""

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np

import openmc
import openmc.examples


def digest_tracks(path):
    tracks = openmc.Tracks(path)
    sha = hashlib.sha256()
    n_tracks = 0
    n_states = 0
    for track in tracks:
        for ptrack in track:
            sha.update(np.int64(ptrack.particle.pdg_number).tobytes())
            # Hash field by field: the raw bytes of a structured array include
            # padding whose layout depends on the HDF5 build.
            for name in ptrack.states.dtype.names:
                field = ptrack.states[name]
                sha.update(np.ascontiguousarray(field, dtype=field.dtype.newbyteorder('<')).tobytes())
            n_tracks += 1
            n_states += len(ptrack.states)
    return {'sha256': sha.hexdigest(), 'tracks': n_tracks, 'states': n_states}


def digest_statepoint(path):
    sha = hashlib.sha256()
    with openmc.StatePoint(path) as sp:
        sha.update(np.ascontiguousarray(sp.k_generation).tobytes())
        for tally in sp.tallies.values():
            sha.update(np.ascontiguousarray(tally.mean).tobytes())
            sha.update(np.ascontiguousarray(tally.std_dev).tobytes())
        k = [float(sp.k_generation[i]) for i in range(len(sp.k_generation))]
    return {'sha256': sha.hexdigest(), 'k_generation': k}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('output')
    parser.add_argument('--particles', type=int, default=2000)
    parser.add_argument('--batches', type=int, default=12)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()

    model = openmc.examples.pwr_pin_cell()
    model.settings.particles = args.particles
    model.settings.batches = args.batches
    model.settings.inactive = 2
    model.settings.seed = 20261017
    model.settings.max_tracks = 10**9
    model.settings.output = {'summary': False}

    # Track every history of the first active batch to keep the file small
    # while still covering fission, scattering, and boundary crossings.
    model.settings.track = [(3, 1, i + 1) for i in range(200)]

    cell_filter = openmc.CellFilter(list(model.geometry.get_all_cells().values()))
    energy_filter = openmc.EnergyFilter(np.logspace(-3, 7, 51))
    tally = openmc.Tally(name='flux')
    tally.filters = [cell_filter, energy_filter]
    tally.scores = ['flux', 'absorption', 'fission']
    model.tallies = openmc.Tallies([tally])

    workdir = Path('track-digest-run')
    workdir.mkdir(exist_ok=True)
    statepoint = model.run(cwd=workdir, threads=args.threads,
                           output=False)

    result = {
        'platform': platform.platform(),
        'machine': platform.machine(),
        'openmc_version': openmc.__version__,
        'tracks': digest_tracks(workdir / 'tracks.h5'),
        'statepoint': digest_statepoint(statepoint),
    }
    with open(args.output, 'w') as fh:
        json.dump(result, fh, indent=2)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
