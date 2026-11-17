"""Watt read coverage using an explicitly synthetic HDF5 law on a real file.

The configured ENDF/B-VIII.0 library has no native WattEnergy entries. This test
uses its H1 file only as a loadable container, and never runs transport.
"""
import shutil

import h5py
import numpy as np

import openmc
import openmc.lib


def test_watt_native_read(tmp_path, mpi_intracomm):
    installed = openmc.data.DataLibrary.from_xml()
    path = tmp_path / 'H1.h5'
    shutil.copyfile(installed.get_by_material('H1')['path'], path)
    a = openmc.data.Tabulated1D([1.0, 1.e6, 2.e7], [5.e5, 8.e5, 1.e6],
                               [2, 3], [1, 2])
    b = openmc.data.Tabulated1D([1.0, 2.e7], [1.e-6, 2.e-6], [2], [2])
    watt = openmc.data.WattEnergy(a, b, -1.e6)
    with h5py.File(path, 'r+') as f:
        d = f['H1/reactions/reaction_102/product_0/distribution_0']
        assert d.attrs['type'] == b'uncorrelated'
        del d['energy']
        watt.to_hdf5(d.create_group('energy'))
    library = openmc.data.DataLibrary()
    library.register_file(path)
    xs = tmp_path / 'cross_sections.xml'
    library.export_to_xml(xs)
    material = openmc.Material()
    material.add_nuclide('H1', 1.0)
    material.set_density('g/cm3', 1.0)
    sphere = openmc.Sphere(r=1, boundary_type='vacuum')
    cell = openmc.Cell(fill=material, region=-sphere)
    model = openmc.Model(geometry=openmc.Geometry([cell]),
                         materials=openmc.Materials([material]))
    model.materials.cross_sections = xs
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    model.settings.batches = 1
    model.settings.particles = 10
    model.settings.output = {'path': str(tmp_path), 'summary': False}
    try:
        model.init_lib(directory=str(tmp_path), output=False, intracomm=mpi_intracomm)
        descriptor = openmc.lib.get_reaction_products('H1', 102)[0]['distributions'][0]
        assert descriptor['energy_type'] == 'watt'
        assert descriptor['read_supported']
        assert not descriptor['energy_mutation_supported']
        actual = openmc.lib.get_emission_distribution('H1', 102)['energy']
        assert actual['type'] == 'watt'
        assert actual['u'] == -1.e6
        for key, expected in (('a', a), ('b', b)):
            for name in ('x', 'y', 'breakpoints', 'interpolation'):
                np.testing.assert_array_equal(
                    actual[key][name], getattr(expected, name))
    finally:
        model.finalize_lib()
