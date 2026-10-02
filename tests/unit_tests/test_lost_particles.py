from pathlib import Path

import openmc
import pytest

from tests.testing_harness import config


@pytest.fixture
def model():
    mat = openmc.Material()
    mat.add_nuclide('N14', 1.0)
    mat.set_density('g/cm3', 1e-5)

    s1 = openmc.Sphere(r=80.0)
    s2 = openmc.Sphere(r=90.0)
    s3 = openmc.Sphere(r=100.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=mat, region=-s1)
    cell2 = openmc.Cell(fill=mat, region=+s2 & -s3)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 50
    model.settings.max_lost_particles = 1000
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())

    return model


def test_max_write_lost_particles(model: openmc.Model, run_in_tmpdir):
    # Set maximum number of lost particle restart files
    model.settings.max_write_lost_particles = 5

    # Run OpenMC to generate lost particle files. Use one thread so that we know
    # exactly how much will be produced. If running in MPI mode, setup proper
    # keyword arguments for run()
    kwargs = {'openmc_exec': config['exe']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    model.run(threads=1, **kwargs)

    # Make sure number of lost particle files is as expected
    lost_particle_files = list(Path.cwd().glob('particle*.h5'))
    n_procs = int(config['mpi_np']) if config['mpi'] else 1
    assert len(lost_particle_files) == model.settings.max_write_lost_particles * n_procs



def test_split_on_coincident_surfaces(run_in_tmpdir):
    """Particles split by weight windows on a surface that coincides with
    another surface must not be lost (GitHub issue #4140)"""
    s1 = openmc.XPlane(-57452.33336021505)
    s2 = openmc.XPlane(43403.32187479845)
    s3 = openmc.XPlane(62125.04976135607)
    s4 = openmc.XPlane(62125.04976135607)
    s5 = openmc.XPlane(100000.0, boundary_type='vacuum')
    inner = openmc.Universe(cells=[
        openmc.Cell(region=-s1 | -s2 | -s3),
        openmc.Cell(region=+s3 & -s5),
    ])
    root = openmc.Cell(fill=inner, region=-s1 | -s2 | -s3 | +s4)

    model = openmc.Model()
    model.geometry = openmc.Geometry([root])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1
    model.settings.batches = 1
    model.settings.max_lost_particles = 100
    model.settings.rel_max_lost_particles = 0.999999
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((-85626.45049221347, 0.0, 0.0)),
        angle=openmc.stats.Monodirectional(
            (0.9182609440159194, 0.39597580569397484, 0.0)),
        energy=openmc.stats.delta_function(1.0e6),
    )

    # The particle is split into ten when crossing s3 into the second mesh bin
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-1.0e5, -1.0e5, -1.0)
    mesh.upper_right = (2.2e5, 1.0e5, 1.0)
    mesh.dimension = (2, 1, 1)
    model.settings.weight_windows = openmc.WeightWindows(
        mesh, lower_ww_bounds=[0.5, 0.01], upper_ww_bounds=[5.0, 0.1],
        energy_bounds=[0.0, 2.0e7])
    model.settings.weight_window_checkpoints = {
        'collision': False, 'surface': True}

    sp_path = model.run()

    # All weight should leak out of the vacuum boundary
    with openmc.StatePoint(sp_path) as sp:
        leakage = sp.global_tallies[sp.global_tallies['name'] == b'leakage']
        assert leakage['mean'][0] == pytest.approx(1.0)
