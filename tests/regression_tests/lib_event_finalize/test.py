"""Event-based particles must not retain caches across library finalization."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import openmc
import openmc.lib


def _make_model(multigroup, directory):
    material = openmc.Material()
    if multigroup:
        groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
        xs = openmc.XSdata('fuel', groups)
        xs.order = 0
        xs.set_total([1.0])
        xs.set_absorption([0.5])
        xs.set_fission([0.15])
        xs.set_nu_fission([0.375])
        xs.set_scatter_matrix([[[0.5]]])
        xs.set_chi([1.0])
        library = openmc.MGXSLibrary(groups)
        library.add_xsdata(xs)
        library.export_to_hdf5(directory / 'mg.h5')
        material.set_density('macro', 1.0)
        material.add_macroscopic('fuel')
    else:
        material.add_nuclide('U235', 2.2e-4)
        material.add_nuclide('H1', 6.6e-2)
        material.set_density('sum')

    materials = openmc.Materials([material])
    if multigroup:
        materials.cross_sections = str(directory / 'mg.h5')
    box = openmc.model.RectangularParallelepiped(
        -5, 5, -5, 5, -5, 5, boundary_type='reflective')
    geometry = openmc.Geometry([openmc.Cell(fill=material, region=-box)])
    model = openmc.Model(geometry=geometry, materials=materials)
    model.settings.energy_mode = (
        'multi-group' if multigroup else 'continuous-energy')
    model.settings.particles = 500
    model.settings.batches = 4
    model.settings.inactive = 1
    model.settings.seed = 1
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-5, -5, -5), (5, 5, 5)))
    return model


def _run_sequence(with_multigroup):
    base = Path.cwd()
    sequence = (True, False) if with_multigroup else (False,)
    for index, multigroup in enumerate(sequence):
        directory = base / str(index)
        directory.mkdir()
        model = _make_model(multigroup, directory)
        os.chdir(directory)
        model.init_lib(output=False, event_based=True)
        try:
            openmc.lib.run(output=False)
            keff = openmc.lib.keff()
        finally:
            model.finalize_lib()
            os.chdir(base)
    (base / 'keff.json').write_text(json.dumps(keff))


def test_lib_event_finalize(tmp_path):
    # Run in subprocesses so a stale particle cache produces an assertion
    # failure rather than terminating the entire pytest session (issue #4200).
    results = []
    for sequence in ('ce', 'mg-ce'):
        directory = tmp_path / sequence
        directory.mkdir()
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), sequence],
            cwd=directory, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, (
            f'{sequence} exited with {result.returncode}\n'
            f'{result.stdout}\n{result.stderr}')
        results.append(json.loads((directory / 'keff.json').read_text()))

    # Both the mean and uncertainty must match a fresh CE run. The same seed
    # and particle count also force reuse of every old MG particle if cleanup
    # is omitted, exposing its empty continuous-energy cross section cache.
    np.testing.assert_allclose(results[1], results[0], rtol=1e-12, atol=0.0)


if __name__ == '__main__':
    _run_sequence(sys.argv[1] == 'mg-ce')
