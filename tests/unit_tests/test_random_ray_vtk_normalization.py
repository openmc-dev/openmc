import math
import struct

import openmc
from openmc.examples import random_ray_three_region_cube
from tests.regression_tests import config


def _read_vtk_flux(path, group, count):
    with open(path, 'rb') as fh:
        data = fh.read()

    marker = (f'SCALARS flux_group_{group} float\n'
              'LOOKUP_TABLE default\n').encode()
    start = data.index(marker) + len(marker)
    end = start + count * struct.calcsize('>f')
    return struct.unpack(f'>{count}f', data[start:end])


def test_random_ray_vtk_normalization(run_in_tmpdir):
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    model.settings.random_ray['source_shape'] = 'linear'
    model.settings.inactive = 5
    model.settings.batches = 10

    plot = openmc.VoxelPlot()
    plot.filename = 'random_ray_vtk_normalization.h5'
    plot.origin = (2.5, 2.5, 2.5)
    plot.width = (5.0, 5.0, 5.0)
    plot.pixels = (2, 2, 2)
    model.plots = [plot]

    model.run(openmc_exec=config['exe'])

    with openmc.StatePoint('statepoint.10.h5') as statepoint:
        source_tally = statepoint.get_tally(name='Source Tally')
        tally_flux = source_tally.mean.ravel()[0]

    vtk_flux = _read_vtk_flux(
        'random_ray_vtk_normalization.vtk', group=0, count=8)
    plot_flux = sum(vtk_flux) / len(vtk_flux)

    assert math.isclose(plot_flux, tally_flux, rel_tol=0.15)
