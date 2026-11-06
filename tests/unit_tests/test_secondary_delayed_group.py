"""Test that split and (n,xn) neutrons keep the delayed group of the neutron
that produced them in a transport run."""

import h5py
import numpy as np
import openmc
import pytest

from tests.regression_tests import config

DELAYED_GROUP = 3
NEUTRON = openmc.ParticleType.NEUTRON.pdg_number
PHOTON = openmc.ParticleType.PHOTON.pdg_number


def count_delayed_groups(sites):
    """Number of sites of each delayed group"""
    groups, counts = np.unique(sites['delayed_group'], return_counts=True)
    return dict(zip(groups.tolist(), counts.tolist()))


@pytest.mark.parametrize('shared_secondary_bank', [False, True],
                         ids=['local_bank', 'shared_bank'])
@pytest.mark.parametrize('energy, weight_windows_on', [
    (14.0e6, False),
    (1.0e6, True),
], ids=['n2n', 'split'])
def test_secondary_neutrons_keep_delayed_group(
        energy, weight_windows_on, shared_secondary_bank, run_in_tmpdir):
    """Start every history from a source file site of delayed group 3 in a
    beryllium sphere with a thin iron shell, and check the delayed group of
    the particles that leave it, which the surface source file records.

    At 14 MeV without weight windows, the only neutrons banked during
    transport are the extra neutrons of (n,2n) reactions, which have integral
    yields in Be9 and Fe56. At 1 MeV, below the (n,2n) thresholds, weight
    windows split the neutrons that collide outside the central mesh cell.
    Nothing fissions, so every neutron continues the history of a source
    neutron and must keep its delayed group, while the photons, mostly from
    inelastic scattering in iron, must have none.
    """
    model = openmc.Model()

    beryllium = openmc.Material()
    beryllium.add_nuclide('Be9', 1.0)
    beryllium.set_density('g/cm3', 1.85)
    iron = openmc.Material()
    iron.add_nuclide('Fe56', 1.0)
    iron.set_density('g/cm3', 7.87)
    radius = 10.0
    inner = openmc.Sphere(r=radius - 1.0)
    outer = openmc.Sphere(r=radius, boundary_type='vacuum')
    model.geometry = openmc.Geometry([
        openmc.Cell(region=-inner, fill=beryllium),
        openmc.Cell(region=+inner & -outer, fill=iron),
    ])

    openmc.write_source_file([openmc.SourceParticle(
        E=energy, delayed_group=DELAYED_GROUP)], 'source.h5')

    n_histories = 1000
    max_sites = 100 * n_histories
    model.settings.run_mode = 'fixed source'
    model.settings.particles = n_histories // 2
    model.settings.batches = 2
    model.settings.photon_transport = True
    model.settings.source = openmc.FileSource('source.h5')
    model.settings.shared_secondary_bank = shared_secondary_bank
    model.settings.surf_source_write = {
        'surface_ids': [outer.id], 'max_particles': max_sites}

    if weight_windows_on:
        # A neutron that first collides in the central cell, which has no
        # window, keeps the nominal windows of the other cells, whose upper
        # bound of 0.1 then splits it into 10 neutrons
        mesh = openmc.RegularMesh()
        mesh.lower_left = (-radius, -radius, -radius)
        mesh.upper_right = (radius, radius, radius)
        mesh.dimension = (3, 3, 3)
        has_window = np.ones(mesh.dimension, dtype=bool)
        has_window[1, 1, 1] = False
        model.settings.weight_windows = [openmc.WeightWindows(
            mesh,
            lower_ww_bounds=np.where(has_window, 0.05, -1.0),
            upper_ww_bounds=np.where(has_window, 0.1, -1.0),
            energy_bounds=[0.0, 2.0e7],
            particle_type='neutron',
        )]
        model.settings.weight_windows_on = True

    kwargs = {'openmc_exec': config['exe'], 'event_based': config['event']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    model.run(**kwargs)

    with h5py.File('surface_source.h5', 'r') as fh:
        sites = fh['source_bank'][...]
    assert len(sites) < max_sites
    neutrons = sites[sites['particle'] == NEUTRON]
    photons = sites[sites['particle'] == PHOTON]

    # A neutron leaves the sphere at most once, so more neutrons than histories
    # means that neutrons banked during transport left it too
    assert len(neutrons) > n_histories
    assert len(photons) > 0

    assert count_delayed_groups(neutrons) == {DELAYED_GROUP: len(neutrons)}
    assert count_delayed_groups(photons) == {0: len(photons)}
