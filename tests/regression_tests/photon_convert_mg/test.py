import numpy as np
import openmc


def test_photon_convert_mg(tmp_path):
    material = openmc.Material(name='lead')
    material.set_density('g/cm3', 11.35)
    material.add_element('Pb', 1.0)

    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=material, region=-sphere)
    model = openmc.Model(
        geometry=openmc.Geometry([cell]),
        materials=openmc.Materials([material]))

    model.settings.run_mode = 'fixed source'
    model.settings.photon_transport = True
    model.settings.atomic_relaxation = True
    model.settings.electron_treatment = 'ttb'
    model.settings.cutoff = {'energy_photon': 1000.0}
    model.settings.source = openmc.IndependentSource(
        particle='photon',
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))

    mgxs_path = tmp_path / 'mgxs.h5'
    model.convert_to_multigroup(
        groups=[1.0e3, 1.0e5, 5.0e5, 1.1e6],
        mgxs_path=mgxs_path,
        particles=2000,
        batches=5,
        inactive=0)

    library = openmc.MGXSLibrary.from_hdf5(mgxs_path)
    assert library.particle_type == openmc.ParticleType.PHOTON
    assert model.settings.energy_mode == 'multi-group'
    assert model.materials[0]._macroscopic

    xsdata = library.xsdatas[0]
    assert np.all(xsdata.total[0] > 0.0)
    assert np.all(xsdata.absorption[0] >= 0.0)
    assert np.any(xsdata.scatter_matrix[0][:, :, 0] > 0.0)
    assert xsdata.multiplicity_matrix[0] is None
