import os
from pathlib import Path

import numpy as np
import openmc

from tests.regression_tests import config


def _photon_model():
    material = openmc.Material(name='lead')
    material.set_density('g/cm3', 11.35)
    material.add_element('Pb', 1.0)

    left = openmc.XPlane(x0=0.0, boundary_type='vacuum')
    right = openmc.XPlane(x0=2.0, boundary_type='vacuum')
    bottom = openmc.YPlane(y0=-0.5, boundary_type='reflective')
    top = openmc.YPlane(y0=0.5, boundary_type='reflective')
    back = openmc.ZPlane(z0=-0.5, boundary_type='reflective')
    front = openmc.ZPlane(z0=0.5, boundary_type='reflective')
    region = +left & -right & +bottom & -top & +back & -front

    model = openmc.Model()
    model.geometry = openmc.Geometry([openmc.Cell(fill=material, region=region)])
    model.materials = openmc.Materials([material])

    mesh = openmc.RegularMesh()
    mesh.dimension = (8, 1, 1)
    mesh.lower_left = (0.0, -0.5, -0.5)
    mesh.upper_right = (2.0, 0.5, 0.5)

    group_edges = [1.0e3, 1.0e5, 5.0e5, 1.1e6]
    tally = openmc.Tally(name='spatial photon flux')
    tally.filters = [
        openmc.MeshFilter(mesh),
        openmc.EnergyFilter(group_edges),
        openmc.ParticleFilter('photon'),
    ]
    tally.scores = ['flux']
    model.tallies.append(tally)

    model.settings.run_mode = 'fixed source'
    model.settings.particles = 10_000
    model.settings.batches = 10
    model.settings.photon_transport = True
    model.settings.atomic_relaxation = True
    model.settings.electron_treatment = 'ttb'
    model.settings.cutoff = {'energy_photon': 1.0e3}
    model.settings.source = openmc.IndependentSource(
        particle='photon',
        space=openmc.stats.Point((0.125, 0.0, 0.0)),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))

    return model, mesh, group_edges


def _spatial_flux(statepoint_path, num_spatial_bins, num_groups):
    with openmc.StatePoint(statepoint_path) as statepoint:
        tally = statepoint.get_tally(name='spatial photon flux')
        flux = tally.mean[:, 0, 0].reshape(num_spatial_bins, num_groups)
    return flux


def test_random_ray_photon_end_to_end(tmp_path, monkeypatch):
    model, mesh, group_edges = _photon_model()
    executable = Path(config['exe'])
    if executable.parent != Path('.'):
        executable = executable.resolve()
        monkeypatch.setenv(
            'PATH', f'{executable.parent}{os.pathsep}{os.environ["PATH"]}')

    ce_directory = tmp_path / 'continuous-energy'
    ce_directory.mkdir()
    ce_statepoint = model.run(cwd=ce_directory, openmc_exec=executable)
    ce_flux = _spatial_flux(
        ce_statepoint, mesh.n_elements, len(group_edges) - 1)

    mgxs_path = tmp_path / 'mgxs.h5'
    model.convert_to_multigroup(
        groups=group_edges,
        particle_type='photon',
        correction=None,
        particles=10_000,
        batches=10,
        overwrite_mgxs_library=True,
        mgxs_path=mgxs_path,
    )
    model.convert_to_random_ray()
    assert model.settings.random_ray['ray_source'].particle \
        == openmc.ParticleType.PHOTON

    model.settings.random_ray['source_region_meshes'] = [
        (mesh, [model.geometry.root_universe])]
    model.settings.random_ray['distance_inactive'] = 50.0
    model.settings.random_ray['distance_active'] = 500.0
    model.settings.inactive = 20
    model.settings.batches = 40
    model.settings.particles = 400

    rr_directory = tmp_path / 'random-ray'
    rr_directory.mkdir()
    rr_statepoint = model.run(cwd=rr_directory, openmc_exec=executable)
    rr_flux = _spatial_flux(
        rr_statepoint, mesh.n_elements, len(group_edges) - 1)

    ce_profile = ce_flux.sum(axis=1) / ce_flux.sum()
    rr_profile = rr_flux.sum(axis=1) / rr_flux.sum()

    cosine_similarity = np.dot(ce_profile, rr_profile) / (
        np.linalg.norm(ce_profile) * np.linalg.norm(rr_profile))
    mean_bin_error = np.mean(np.abs(ce_profile - rr_profile))
    ce_group_fractions = ce_flux.sum(axis=0) / ce_flux.sum()
    rr_group_fractions = rr_flux.sum(axis=0) / rr_flux.sum()
    integral_ratio = rr_flux.sum() / ce_flux.sum()

    assert cosine_similarity > 0.97
    assert mean_bin_error < 0.02
    assert np.max(np.abs(ce_group_fractions - rr_group_fractions)) < 0.03
    assert rr_group_fractions[:-1].sum() > 0.05
    assert 0.8 < integral_ratio < 1.2
