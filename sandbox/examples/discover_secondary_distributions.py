#!/usr/bin/env python3
"""Discover native reaction secondary laws without running or modifying transport.

Run with OPENMC_CROSS_SECTIONS pointing to an installed evaluated library:
    python -m sandbox.examples.discover_secondary_distributions
    python -m sandbox.examples.discover_secondary_distributions Li6 Be9 Fe56

Li6 spectra are saved as PNG files in li6-spectra. Tabulated laws use the
nearest incident-energy table; line spectra are conditional probabilities,
not reaction yields. Elastic spectra assume a stationary target.
"""
import argparse
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import matplotlib.pyplot as plt
import numpy as np

import openmc
import openmc.lib


def plot_li6_spectra(nuclide_name, products, library, directory, incident_energy):
    """Plot each available Li6 law without modifying native distributions."""
    directory.mkdir(parents=True, exist_ok=True)
    with h5py.File(library.get_by_material(nuclide_name)['path'], 'r') as data:
        mass = float(data[nuclide_name].attrs['atomic_weight_ratio'])
    for product in products:
        for descriptor in product['distributions']:
            mt = product['mt']
            pi = product['product_index']
            di = descriptor['distribution_index']
            law = openmc.lib.get_emission_distribution(nuclide_name, mt, pi, di)
            energy = law.get('energy') if law['type'] == 'uncorrelated' else law
            energy_in = incident_energy
            fig, ax = plt.subplots()
            frame = product['reference_frame']
            print("product ", product.keys())
            kind = energy['type'] if energy is not None else None
            if kind in ('continuous-tabular', 'correlated', 'kalbach-mann'):
                i = int(np.argmin(abs(energy['incident_energy'] - energy_in)))
                energy_in = energy['incident_energy'][i]
                table = energy['tables'][i]
                print("energy keys ", energy.keys())
                start = table['n_discrete']
                print("table keys ", table.keys())
                ax.plot(table['energy'][start:] / 1.e6,
                        table['pdf'][start:] * 1.e6,
                        drawstyle=('steps-post' if table['interpolation'] == 1
                                   else 'default'))
                ax.set_ylabel('Probability density [1/MeV]')
                if start:
                    lines = ax.twinx()
                    lines.vlines(table['energy'][:start] / 1.e6, 0,
                                 table['pdf'][:start], color='tab:orange')
                    lines.set_ylabel('Discrete probability')
            elif kind in ('discrete-photon', 'level-inelastic'):
                if kind == 'discrete-photon':
                    outgoing = energy['energy']
                    if energy['primary_flag'] == 2:
                        a = energy['atomic_weight_ratio']
                        outgoing += a / (a + 1.) * energy_in
                    frame = 'laboratory'
                else:
                    if energy_in <= energy['threshold']:
                        plt.close(fig)
                        print(f'{nuclide_name} MT={mt}: below inelastic threshold; skipped')
                        continue
                    outgoing = energy['mass_ratio'] * (
                        energy_in - energy['threshold'])
                ax.vlines(outgoing / 1.e6, 0, 1)
                ax.set_ylim(0, 1.1)
                ax.set_ylabel('Discrete probability')
            elif mt == 2 and energy is None:
                angle = law['angle']
                if angle['type'] == 'isotropic':
                    mu, pdf, interpolation = np.array([-1., 1.]), np.array([.5, .5]), 1
                elif angle['type'] == 'tabular-angle':
                    i = int(np.argmin(abs(angle['incident_energy'] - energy_in)))
                    energy_in = angle['incident_energy'][i]
                    table = angle['tables'][i]
                    mu, pdf = table['x'], table['pdf']
                    interpolation = table['interpolation']
                else:
                    plt.close(fig)
                    print(f'{nuclide_name} MT={mt}: unsupported angular law; skipped')
                    continue
                jacobian = energy_in * 2. * mass / (mass + 1.)**2
                outgoing = energy_in * (mass**2 + 1. + 2. * mass * mu) / (mass + 1.)**2
                ax.plot(outgoing / 1.e6, pdf / jacobian * 1.e6,
                        drawstyle='steps-post' if interpolation == 1 else 'default')
                frame = 'laboratory (stationary target)'
                ax.set_ylabel('Probability density [1/MeV]')
            else:
                plt.close(fig)
                print(f'{nuclide_name} MT={mt}: unsupported spectrum law {kind}; skipped')
                continue
            ax.set_xlabel(f'Outgoing {product["particle"]} energy [MeV] ({frame})')
            ax.set_title(f'{nuclide_name} MT={mt}, product={pi}, distribution={di}\n'
                         f'Incident energy {energy_in:g} eV')
            fig.tight_layout()
            path = directory / f'{nuclide_name}_mt{mt}_product{pi}_distribution{di}.png'
            fig.savefig(path)
            plt.close(fig)
            print(f'Saved {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('nuclides', nargs='*', default=[
        'H1', 'H2', 'He3', 'He4', 'Li6', 'Li7', 'Be9', 'C12', 'O16',
        'Al27', 'Fe54', 'Fe56', 'Fe57', 'Fe58', 'U235', 'U238', 'Ra223'])
    parser.add_argument('--output', type=Path, default=Path('li6-spectra'))
    parser.add_argument('--incident-energy', type=float, default=2.8e6 +0*14.e6,
                        help='incident neutron energy in eV (default: 14 MeV)')
    args = parser.parse_args()
    if not np.isfinite(args.incident_energy) or args.incident_energy <= 0:
        parser.error('--incident-energy must be finite and positive')
    library = openmc.data.DataLibrary.from_xml()
    material = openmc.Material()
    available = []
    for name in args.nuclides:
        if library.get_by_material(name) is None:
            print(f'{name}: absent from configured library')
            continue
        material.add_nuclide(name, 1.0)
        available.append(name)
    if not available:
        return
    material.set_density('g/cm3', 1.0)
    cell = openmc.Cell(fill=material,
                       region=-openmc.Sphere(r=1, boundary_type='vacuum'))
    model = openmc.Model(geometry=openmc.Geometry([cell]))
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    model.settings.particles = 10
    model.settings.batches = 1
    with TemporaryDirectory() as directory:
        model.settings.output = {'path': directory, 'summary': False}
        try:
            model.init_lib(directory=directory, output=False)
            print('| Nuclide | Example MT/product/distribution | Particle | Mode | '
                  'Frame | AngleEnergy | Energy law | Occurrences | Read | Write |')
            print('|---|---|---|---|---|---|---|---:|---|---|')
            for name in available:
                groups = defaultdict(list)
                products = openmc.lib.get_reaction_products(name)
                for p in products:
                    for d in p['distributions']:
                        key = (p['particle'], p['emission_mode'], p['reference_frame'],
                               d['angle_energy_type'], d['energy_type'])
                        groups[key].append((p, d))
                for key, entries in groups.items():
                    p, d = entries[0]
                    # Read the complete native payload for every reported example.
                    law = openmc.lib.get_emission_distribution(
                        name, p['mt'], p['product_index'], d['distribution_index'])
                    assert law['type'] == d['angle_energy_type']
                    write = ('energy only (legacy API)'
                             if d['energy_mutation_supported'] else 'no')
                    address = (f"{p['mt']}/{p['product_index']}/"
                               f"{d['distribution_index']}")
                    print(f"| {name} | {address} | " + ' | '.join(
                        str(x) if x is not None else '—' for x in key)
                        + f" | {len(entries)} | yes | {write} |")
                if name == 'C12':
                    nuclide_name = name
                    plot_li6_spectra(nuclide_name, products, library, args.output,
                                     args.incident_energy)
        finally:
            model.finalize_lib()

if __name__ == '__main__':
    main()
