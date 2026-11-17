#!/usr/bin/env python3
"""Plot NUCLIDE_OF_INTEREST angular PDFs copied from OpenMC's loaded native objects.

Run with OPENMC_CROSS_SECTIONS configured and the runtime library built:
    python -m sandbox.examples.plot_angular_distributions

No transport is run. All angular data come from openmc.lib snapshots, without
reading evaluated-data files in Python. The nearest incident-energy table is
shown, rather than mixing adjacent tables at the requested energy. Reference
frames are native product frames; no laboratory-frame conversion is applied.
"""
import argparse
from pathlib import Path
from pprint import pformat
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np

import openmc
import openmc.lib


INTERPOLATION = {1: 'histogram', 2: 'linear-linear'}
NUCLIDE_OF_INTEREST = 'C12'
ENERGY_OF_INTEREST = 2e6 + 0*14.e6


def table_integral(table):
    """Integrate a native continuous angular table using its interpolation."""
    # AngleDistribution stores Tabular objects with no discrete angular masses.
    # n_discrete in emission snapshots belongs to outgoing ENERGY tables only.
    x, pdf = table['x'], table['pdf']
    widths = np.diff(x)
    if table['interpolation'] == 1:
        return float(np.sum(widths * pdf[:-1]))
    if table['interpolation'] == 2:
        return float(np.sum(widths * (pdf[:-1] + pdf[1:]) / 2))
    raise ValueError(f"Unsupported native interpolation {table['interpolation']}")


def compact(value):
    """Abbreviate arrays and repeated tables while retaining dictionary keys."""
    if isinstance(value, np.ndarray):
        return f'ndarray(shape={value.shape}, dtype={value.dtype}, {value!r})'
    if isinstance(value, dict):
        return {key: compact(item) for key, item in value.items()}
    if isinstance(value, list) and len(value) > 2:
        return [compact(value[0]), f'... {len(value) - 2} tables omitted ...',
                compact(value[-1])]
    if isinstance(value, list):
        return [compact(item) for item in value]
    return value


def theta_curve(table):
    """Transform each native mu interval into degrees without smoothing jumps."""
    x, pdf = table['x'], table['pdf']
    theta_parts, density_parts = [], []
    # Reverse interval order so theta increases from forward to backward.
    for i in range(len(x) - 2, -1, -1):
        theta = np.linspace(np.arccos(x[i + 1]), np.arccos(x[i]), 65)
        mu = np.cos(theta)
        if table['interpolation'] == 1:
            density = np.full(theta.shape, pdf[i])
        elif table['interpolation'] == 2:
            fraction = (mu - x[i]) / (x[i + 1] - x[i])
            density = pdf[i] + fraction * (pdf[i + 1] - pdf[i])
        else:
            raise ValueError('Unsupported native angular interpolation')
        # A histogram in mu becomes a sine curve within each theta interval,
        # not a histogram in theta. NaNs preserve the native bin discontinuities.
        theta_parts.extend((np.degrees(theta), [np.nan]))
        density_parts.extend((density * np.sin(theta) * np.pi / 180, [np.nan]))
    return np.concatenate(theta_parts), np.concatenate(density_parts)


def plot_law(product, descriptor, law, requested_energy, directory):
    """Plot the independent angular component and report its native metadata."""
    mt, pi = product['mt'], product['product_index']
    di = descriptor['distribution_index']
    address = f'{NUCLIDE_OF_INTEREST} MT={mt} product={pi} distribution={di}'
    if law['type'] in ('correlated', 'kalbach-mann', 'nbody'):
        print(f'{address}: angular law {law["type"]} is conditional/correlated; '
              'no independent p(mu | E_in) plotted (no marginalisation).')
        return False
    if law['type'] != 'uncorrelated':
        print(f'{address}: unsupported law {law["type"]}; skipped')
        return False
    angle = law['angle']
    selected_energy = requested_energy
    if angle['type'] == 'tabular-angle':
        index = int(np.argmin(abs(angle['incident_energy'] - requested_energy)))
        selected_energy = float(angle['incident_energy'][index])
        table = angle['tables'][index]
        print(f'{address}: selected angular table {index} of '
              f'{len(angle["tables"])}; incident sampling interpolation = '
              f'{angle["sampling_interpolation"]}')
    elif angle['type'] == 'isotropic':
        # Exact analytic rendering of the native isotropic tag (no stored grid).
        table = {'x': np.array([-1., 1.]), 'pdf': np.array([.5, .5]),
                 'interpolation': 1}
        print(f'{address}: isotropic at all energies; no incident grid')
    else:
        print(f'{address}: unsupported angular representation {angle["type"]}')
        return False
    interpolation = table['interpolation']
    if interpolation not in INTERPOLATION:
        print(f'{address}: unsupported native interpolation {interpolation}; skipped')
        return False
    mass = table_integral(table)
    if mt == 2:
        with np.printoptions(threshold=8, edgeitems=2):
            print('{NUCLIDE_OF_INTEREST} MT=2 selected native angular table:')
            print(pformat(table, sort_dicts=False))
    print(f'{address}\n'
          f'  emitted particle = {product["particle"]}\n'
          f'  reference frame = {product["reference_frame"]}\n'
          f'  angular representation = {angle["type"]}\n'
          f'  incident energy = {selected_energy:.6e} eV\n'
          f'  interpolation = {interpolation} ({INTERPOLATION[interpolation]})\n'
          f'  integral p(mu) dmu = {mass:.12f}; discrete angular mass = 0\n'
          f'  normalisation check = '
          f'{"PASS" if np.isclose(mass, 1., rtol=0., atol=1.e-6) else "FAIL"}')
    if 'integral' in table:
        print(f'  native loader integral (before loader normalisation) = '
              f'{table["integral"]:.12f}')
    title = (f'{NUCLIDE_OF_INTEREST} MT={mt}, emitted {product["particle"]} '
             f'(product {pi}, distribution {di})\n'
             f'E_in = {selected_energy:.6e} eV; '
             f'{product["reference_frame"]}')
    stem = (f'{NUCLIDE_OF_INTEREST}_mt2' if (mt, pi, di) == (2, 0, 0) else
            f'{NUCLIDE_OF_INTEREST}_mt{mt}_product{pi}_distribution{di}')
    for coordinate in ('mu', 'theta'):
        fig, ax = plt.subplots(figsize=(8, 5))
        if coordinate == 'mu':
            ax.plot(table['x'], table['pdf'], linewidth=1.8,
                    drawstyle='steps-post' if interpolation == 1 else 'default')
            ax.set_xlim(-1, 1)
            ax.set_xlabel('Scattering cosine, μ = cos(θ)')
            ax.set_ylabel('Probability density, p(μ)')
        else:
            theta, density = theta_curve(table)
            ax.plot(theta, density, linewidth=1.8)
            ax.set_xlim(0, 180)
            ax.set_xlabel('Scattering angle, θ [degrees] (0° = forward)')
            ax.set_ylabel('Probability density [1/degree]')
        ax.set_ylim(bottom=0)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=.25)
        fig.tight_layout()
        filename = directory / f'{stem}_{coordinate}.png'
        fig.savefig(filename, dpi=180)
        plt.close(fig)
        print(f'  output PNG = {filename}')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--incident-energy', type=float, default=ENERGY_OF_INTEREST,
                        help='requested incident energy in eV (default: ENERGY_OF_INTEREST)')
    parser.add_argument('--output', type=Path,
                        default=Path('angular-distribution-plots'))
    args = parser.parse_args()
    if not np.isfinite(args.incident_energy) or args.incident_energy <= 0:
        parser.error('--incident-energy must be finite and positive')
    args.output.mkdir(parents=True, exist_ok=True)
    material = openmc.Material()
    material.add_nuclide(NUCLIDE_OF_INTEREST, 1.)
    material.set_density('g/cm3', 1.)
    sphere = openmc.Sphere(r=1., boundary_type='vacuum')
    model = openmc.Model(geometry=openmc.Geometry([
        openmc.Cell(fill=material, region=-sphere)]))
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    model.settings.particles = 10
    model.settings.batches = 1
    with TemporaryDirectory(prefix='openmc-angular-') as directory:
        model.settings.output = {'path': directory, 'summary': False}
        try:
            model.init_lib(directory=directory, output=False)
            products = openmc.lib.get_reaction_products(NUCLIDE_OF_INTEREST)
            print(f'Discovered {NUCLIDE_OF_INTEREST} products (zero-based product/distribution indices):')
            for product in products:
                laws = ', '.join(
                    f'{d["distribution_index"]}:{d["angle_energy_type"]}'
                    for d in product['distributions']) or '(none)'
                print(f'  MT={product["mt"]} product={product["product_index"]} '
                      f'particle={product["particle"]} '
                      f'frame={product["reference_frame"]} distributions={laws}')
            elastic_plotted = False
            for product in sorted(products, key=lambda p: (p['mt'] != 2, p['mt'],
                                                          p['product_index'])):
                if product['particle'] != 'neutron':
                    print(f'{NUCLIDE_OF_INTEREST} MT={product["mt"]} '
                          f'product={product["product_index"]}: '
                          f'{product["particle"]} product; neutron example skips it')
                    continue
                for descriptor in product['distributions']:
                    mt, pi = product['mt'], product['product_index']
                    di = descriptor['distribution_index']
                    if not descriptor['read_supported']:
                        print(f'{NUCLIDE_OF_INTEREST} MT={mt} product={pi} distribution={di}: '
                              f'{descriptor["read_reason"]}')
                        continue
                    law = openmc.lib.get_emission_distribution(NUCLIDE_OF_INTEREST, mt, pi, di)
                    if mt == 2:
                        with np.printoptions(threshold=8, edgeitems=2):
                            print(f'{NUCLIDE_OF_INTEREST} MT=2 get_emission_distribution() (abbreviated):')
                            print(pformat(compact(law), sort_dicts=False))
                        print('Angular table schema: '
                              '{type, x, pdf, cdf, interpolation, integral}; '
                              'x/pdf/cdf are independent float64 array copies. '
                              'No discrete angular component in native Tabular.')
                    plotted = plot_law(product, descriptor, law,
                                       args.incident_energy, args.output)
                    elastic_plotted |= mt == 2 and plotted
            if not elastic_plotted:
                raise RuntimeError(f'No {NUCLIDE_OF_INTEREST} elastic angular distribution was plotted')
        finally:
            model.finalize_lib()


if __name__ == '__main__':
    main()
