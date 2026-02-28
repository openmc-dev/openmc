import glob
import os

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness

def delta_tracking_model(nuclide):
    model = openmc.model.Model()

    model.materials = openmc.Materials()

    # single material
    mat = openmc.Material()
    mat.add_nuclide(nuclide, 1.0)
    mat.set_density('g/cc', 19.1)
    mat.depletable = True

    model.materials.append(mat)

    # geometry (infinite cell medium)
    cell = openmc.Cell()
    cell.fill = mat
    cell.temperature = 300

    prism = openmc.model.rectangular_prism(1000.00, 1000.00, boundary_type='vacuum')

    cell.region = prism

    model.geometry = openmc.Geometry([cell,])

    settings = model.settings
    settings.run_mode = 'fixed source'
    settings.particles = 100
    settings.batches = 10
    settings.inactive = 1
    settings.delta_tracking = True

    return model

nuclides = ('Fe56', 'U238', 'H1', 'O16', 'Zr91', 'Zr90', 'Zr92',  'Zr94',  'Zr96', 'B10', 'B11')
@pytest.mark.parametrize('nuclide', nuclides)
def test_delta_tracking(nuclide):
    model = delta_tracking_model(nuclide)
    harness = PyAPITestHarness('statepoint.10.h5',model=model)
    harness._build_inputs()
    harness._run_openmc()
    harness._cleanup()
    output = glob.glob("*.txt")
    for f in output:
        if os.path.exists(f):
            os.remove(f)
