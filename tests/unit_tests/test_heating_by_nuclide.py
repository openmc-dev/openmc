import openmc
import pytest

from tests.testing_harness import config


@pytest.fixture
def model():
    # Create simple sphere model
    model = openmc.Model()
    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.add_nuclide('H1', 1.0)
    mat.set_density('g/cm3', 5.0)
    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 1000
    model.settings.inactive = 0
    model.settings.batches = 5
    model.settings.photon_transport = True

    # Add two tallies, one with heating by nuclide and one with total heating
    particle_filter = openmc.ParticleFilter(['neutron', 'photon'])
    heating_by_nuclide = openmc.Tally()
    heating_by_nuclide.filters = [particle_filter]
    heating_by_nuclide.nuclides = ['U235', 'H1']
    heating_by_nuclide.scores = ['heating']

    heating_total = openmc.Tally()
    heating_total.filters = [particle_filter]
    heating_total.scores = ['heating']
    model.tallies.extend([heating_by_nuclide, heating_total])

    return model


def test_heating_by_nuclide(model, run_in_tmpdir):
    # If running in MPI mode, setup proper keyword arguments for run()
    kwargs = {'openmc_exec': config['exe']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    sp_path = model.run(**kwargs)

    # Get tallies from resulting statepoint
    with openmc.StatePoint(sp_path) as sp:
        heating_by_nuclide = sp.tallies[model.tallies[0].id]
        heating_total = sp.tallies[model.tallies[1].id]

    for particle in heating_by_nuclide.filters[0].bins:
        # Get slice of each tally corresponding to a single particle
        kwargs = {'filters': [openmc.ParticleFilter], 'filter_bins': [(particle,)]}
        particle_slice_by_nuclide = heating_by_nuclide.get_values(**kwargs)
        particle_slice_total = heating_total.get_values(**kwargs)

        # Summing over nuclides should equal total
        assert particle_slice_by_nuclide.sum() == pytest.approx(particle_slice_total.sum())
