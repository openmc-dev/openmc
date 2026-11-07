#!/usr/bin/env python3
"""Transport-level benchmark: system math library versus the CORE-MATH backend.

Runs the same models with two OpenMC executables, one built with the default
system math library and one with OPENMC_USE_COREMATH=ON, and reports the wall
time of each run so the cost of correctly rounded transcendental functions is
measured where it matters, in transport, rather than per call.

Usage:
  coremath-benchmark.py --system /path/to/default/bin/openmc \
      --coremath /path/to/coremath/bin/openmc --repeats 3 --threads 2 \
      --output bench.json
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

import openmc
import openmc.examples


def models():
    pin = openmc.examples.pwr_pin_cell()
    pin.settings.particles = 20000
    pin.settings.batches = 30
    pin.settings.inactive = 5
    pin.settings.seed = 1
    yield 'pwr_pin_cell (CE eigenvalue, 20k x 30)', pin

    asm = openmc.examples.pwr_assembly()
    asm.settings.particles = 10000
    asm.settings.batches = 20
    asm.settings.inactive = 5
    asm.settings.seed = 1
    yield 'pwr_assembly (CE eigenvalue, 10k x 20)', asm

    # Non-multiplying fixed-source problem with coupled neutron-photon
    # transport: a 14.1 MeV neutron point source in a water sphere inside an
    # iron shell, which exercises the photon physics (Compton, pair production,
    # bremsstrahlung) and the scattering kinematics.
    water = openmc.Material(name='water')
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    iron = openmc.Material(name='iron')
    iron.add_element('Fe', 1.0)
    iron.set_density('g/cm3', 7.87)
    inner = openmc.Sphere(r=30.0)
    outer = openmc.Sphere(r=40.0, boundary_type='vacuum')
    cells = [openmc.Cell(fill=water, region=-inner),
             openmc.Cell(fill=iron, region=+inner & -outer)]
    fs = openmc.Model()
    fs.geometry = openmc.Geometry(cells)
    fs.materials = openmc.Materials([water, iron])
    fs.settings.run_mode = 'fixed source'
    fs.settings.particles = 20000
    fs.settings.batches = 10
    fs.settings.photon_transport = True
    fs.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0)),
        energy=openmc.stats.Discrete([14.1e6], [1.0]))
    fs.settings.seed = 1
    tally = openmc.Tally(name='photon flux')
    tally.filters = [openmc.ParticleFilter(['neutron', 'photon']),
                     openmc.EnergyFilter(np.logspace(3, 7.2, 43))]
    tally.scores = ['flux']
    fs.tallies = openmc.Tallies([tally])
    yield 'water sphere, 14 MeV source, neutron-photon (20k x 10)', fs


def run_once(exe, cwd, threads):
    # Each executable lives in its own install prefix; make sure it loads the
    # libopenmc next to it rather than one found elsewhere on the path.
    libdir = str(Path(exe).resolve().parent.parent / 'lib')
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    for var in ('LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH'):
        env[var] = libdir + (':' + env[var] if env.get(var) else '')
    t0 = time.perf_counter()
    proc = subprocess.run([exe, '-s', str(threads)], cwd=cwd, env=env,
                          capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        print(proc.stdout[-3000:])
        print(proc.stderr[-3000:])
        raise RuntimeError(f'{exe} exited with {proc.returncode}')
    return elapsed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--system', required=True)
    parser.add_argument('--coremath', required=True)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--output', default='bench.json')
    args = parser.parse_args()

    results = {}
    for name, model in models():
        with tempfile.TemporaryDirectory() as tmp:
            model.export_to_model_xml(Path(tmp) / 'model.xml')
            times = {'system': [], 'coremath': []}
            for i in range(args.repeats):
                for label, exe in (('system', args.system),
                                   ('coremath', args.coremath)):
                    t = run_once(exe, tmp, args.threads)
                    times[label].append(t)
                    print(f'{name} | {label} | run {i + 1}: {t:.2f} s',
                          flush=True)
            s = statistics.median(times['system'])
            p = statistics.median(times['coremath'])
            results[name] = {
                'system_s': times['system'], 'coremath_s': times['coremath'],
                'system_median_s': s, 'coremath_median_s': p,
                'ratio_coremath_over_system': p / s,
            }
            print(f'{name}: system {s:.2f} s, coremath {p:.2f} s, '
                  f'ratio {p / s:.3f}', flush=True)

    with open(args.output, 'w') as fh:
        json.dump(results, fh, indent=2)
    return 0


if __name__ == '__main__':
    sys.exit(main())
