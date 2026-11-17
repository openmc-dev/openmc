"""Read-only angular API lifecycle, independent of simulation initialization."""
import numpy as np
import pytest

import openmc
import openmc.lib
from openmc import exceptions


def test_angular_lifecycle(tmp_path, mpi_intracomm):
    assert not openmc.lib.is_initialized
    for read in (lambda: openmc.lib.get_angular_distribution(0, 2),
                 lambda: openmc.lib.sample_angular_distribution(0, 2, 1.e6, 10)):
        with pytest.raises(exceptions.OutOfBoundsError):
            read()
    material = openmc.Material()
    material.add_nuclide('Li6', 1.)
    material.set_density('g/cm3', 1.)
    sphere = openmc.Sphere(r=1, boundary_type='vacuum')
    model = openmc.Model(geometry=openmc.Geometry([
        openmc.Cell(fill=material, region=-sphere)]))
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    model.settings.particles = 10
    model.settings.batches = 2
    try:
        model.init_lib(directory=str(tmp_path), output=False, intracomm=mpi_intracomm)
        angle = openmc.lib.get_angular_distribution('Li6', 2)
        saved = angle['tables'][0]['pdf'].copy()
        samples = openmc.lib.sample_angular_distribution('Li6', 2, 1.e6, 10)
    finally:
        model.finalize_lib()
    np.testing.assert_array_equal(angle['tables'][0]['pdf'], saved)
    assert np.all(np.isfinite(samples))
    for read in (lambda: openmc.lib.get_angular_distribution(0, 2),
                 lambda: openmc.lib.sample_angular_distribution(0, 2, 1.e6, 10)):
        with pytest.raises(exceptions.OutOfBoundsError):
            read()
