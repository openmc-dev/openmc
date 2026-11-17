"""Demonstrate in-memory U235 prompt-fission spectrum editing with NNDC data.

Run with OPENMC_CROSS_SECTIONS pointing to NNDC's cross_sections.xml:
    OMP_NUM_THREADS=2 python -m sandbox.examples.runtime_emission_spectrum --output demo

Two runs use identical seeds and an unchanged first batch. After that batch,
only the second run softens the U235 prompt spectrum. Subsequent source-bank
energies demonstrate the effect on transported particles. This is a deliberate
sensitivity perturbation, not an alternative evaluated nuclear-data library.

Also plot Li6 MT=2 using runtime angular data and stationary-target elastic
kinematics at the tabulated incident energy nearest 1 MeV. Elastic scattering
has no separate continuous-tabular energy law for get_emission_spectrum.
"""
from argparse import ArgumentParser
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import h5py
import numpy as np

import openmc
import openmc.lib
from openmc import exceptions as exc


def run_case(directory, modify, particles=2000, batches=5, event_based=False):
    """Run one control/perturbed case and return native tables and bank energies."""
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    openmc.reset_auto_ids()
    model = openmc.examples.pwr_pin_cell()
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.inactive = 0
    model.settings.generations_per_batch = 1
    model.settings.seed = 19
    model.settings.verbosity = 1
    model.settings.output = {'path': str(directory), 'summary': False,
                             'tallies': False}
    model.settings.statepoint = {'batches': []}
    banks = []
    try:
        model.init_lib(directory=str(directory), output=False, event_based=event_based)
        # Load the same nuclides in both cases to preserve identical RNG setup.
        openmc.lib.load_nuclide('Li6')
        if modify:
            plot_li6_elastic(directory.parent / 'li6_mt2_spectrum.png')
        openmc.lib.simulation_init()
        index = openmc.lib.nuclides['U235']._index
        try:
            original = openmc.lib.get_emission_spectrum(index, 18)
        except exc.DataError as err:
            raise RuntimeError(
                'This example requires a continuous-tabular U235 MT=18 prompt '
                'spectrum, as provided by the NNDC HDF5 library.') from err
        # Plot the native spectrum before transport or mutation.
        if modify:
            plot_spectrum(original, directory.parent / 'original.png')
        openmc.lib.next_batch()
        banks.append(openmc.lib.source_bank()['E'].copy())
        if modify:
            changed = deepcopy(original)
            for table in changed['tables']:
                table['pdf'] *= np.exp(-table['energy'] / 13.e5)
            openmc.lib.set_emission_spectrum(index, 18, changed)
        readback = openmc.lib.get_emission_spectrum(index, 18)
        if modify:
            assert any(not np.array_equal(a['pdf'], b['pdf'])
                       for a, b in zip(original['tables'], readback['tables']))
            plot_spectrum(original, directory.parent / 'spectra.png', readback)
        for _ in range(batches - 1):
            openmc.lib.next_batch()
            banks.append(openmc.lib.source_bank()['E'].copy())
        return original, readback, np.asarray(banks)
    finally:
        model.finalize_lib()


def plot_li6_elastic(path, incident_energy=1.e6):
    """Plot Li6 elastic outgoing energies for a stationary target.

    Read native angular tables through the runtime API. The mass ratio, which
    the runtime API does not expose, comes from the configured HDF5 library.
    Use the nearest incident table without interpolating between tables.
    """
    law = openmc.lib.get_emission_distribution('Li6', 2)
    if law['type'] != 'uncorrelated' or law['energy'] is not None:
        raise RuntimeError('Expected an elastic angular-only Li6 MT=2 law')
    angle = law['angle']
    if angle['type'] == 'isotropic':
        energy_in = incident_energy
        mu, pdf, interpolation = np.array([-1., 1.]), np.array([.5, .5]), 1
    elif angle['type'] == 'tabular-angle':
        i = int(np.argmin(np.abs(angle['incident_energy'] - incident_energy)))
        energy_in = angle['incident_energy'][i]
        table = angle['tables'][i]
        mu, pdf = table['x'], table['pdf']
        interpolation = table['interpolation']
    else:
        raise RuntimeError(f'Unsupported elastic angular law: {angle["type"]}')

    library = openmc.data.DataLibrary.from_xml()
    with h5py.File(library.get_by_material('Li6')['path'], 'r') as data:
        mass_ratio = float(data['Li6'].attrs['atomic_weight_ratio'])
    # Elastic transport uses CM angles. Transform the PDF with dE_out/dmu.
    jacobian = energy_in * 2. * mass_ratio / (mass_ratio + 1.)**2
    energy_out = energy_in * (mass_ratio**2 + 1. + 2. * mass_ratio * mu) / (
        mass_ratio + 1.)**2
    fig, ax = plt.subplots()
    ax.plot(energy_out / 1.e6, pdf / jacobian * 1.e6,
            drawstyle='steps-post' if interpolation == 1 else 'default')
    ax.set(xlabel='Outgoing neutron energy [MeV]', ylabel='PDF [1/MeV]',
           title=f'Li6 MT=2 elastic scattering, incident energy {energy_in:g} eV')
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_spectrum(original, path, modified=None):
    """Plot one native outgoing table at the lowest incident energy."""
    fig, ax = plt.subplots()
    for spectrum, label in ((original, 'Original'), (modified, 'Modified')):
        if spectrum is None:
            continue
        table = spectrum['tables'][0]
        start = table['n_discrete']
        ax.plot(table['energy'][start:] / 1.e6,
                table['pdf'][start:] * 1.e6, label=label,
                drawstyle=('steps-post' if table['interpolation'] == 1
                           else 'default'))
    ax.set(xlabel='Outgoing neutron energy [MeV]', ylabel='PDF [1/MeV]',
           xlim=(0, 8), title=(f'U235 prompt fission, incident energy '
                              f'{original["incident_energy"][0]:g} eV'))
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('emission-demo'))
    parser.add_argument('--particles', type=int, default=2000)
    parser.add_argument('--batches', type=int, default=5)
    parser.add_argument('--event', action='store_true')
    args = parser.parse_args()
    if args.batches < 2 or args.particles < 100:
        parser.error('use at least two batches and 100 particles')
    _, _, control = run_case(args.output / 'control', False, args.particles,
                             args.batches, args.event)
    _, _, modified = run_case(args.output / 'modified', True, args.particles,
                              args.batches, args.event)
    np.testing.assert_array_equal(control[0], modified[0])
    if np.array_equal(control[1:], modified[1:]):
        raise RuntimeError('Transport did not respond to the spectrum change')
    np.savez(args.output / 'source_energies.npz', control=control,
             modified=modified)
    means = np.column_stack((np.arange(1, args.batches + 1),
                             control.mean(axis=1) / 1.e6,
                             modified.mean(axis=1) / 1.e6))
    np.savetxt(args.output / 'batch_means.csv', means, delimiter=',',
               header='batch,control_MeV,modified_MeV', comments='')
    fig, ax = plt.subplots()
    bins = np.linspace(0, 8, 81)
    for values, label in ((control, 'Control'), (modified, 'Modified')):
        ax.hist(values[1:].ravel() / 1.e6, bins=bins, density=True,
                histtype='step', label=label)
    ax.set(xlabel='Fission source-bank neutron energy [MeV]',
           ylabel='Probability density [1/MeV]',
           title='Transport after the first batch')
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output / 'transport.png')
    plt.close(fig)
    print(f'Identical first batch confirmed; results saved to {args.output}')
    print(f'Mean source energy after mutation: control '
          f'{control[1:].mean() / 1.e6:.3f} MeV, modified '
          f'{modified[1:].mean() / 1.e6:.3f} MeV')


if __name__ == '__main__':
    main()
