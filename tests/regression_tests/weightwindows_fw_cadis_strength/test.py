import os

import numpy as np
import openmc

from tests.testing_harness import WeightWindowPyAPITestHarness


class MGXSTestHarness(WeightWindowPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_weight_windows_fw_cadis_strength():
    """FW-CADIS weight windows do not depend on the strength of the source.

    A mesh cuts a thin shell holding the source into parts that few rays cross.
    Those that rays only cross in their dead zones are discovered again in the
    adjoint solve, where they must not get the forward source: far from a
    strength of one, it would swamp the adjoint source, which goes as its
    inverse.

    """
    groups = openmc.mgxs.EnergyGroups([1e-5, 20.0e6])
    xsdata = openmc.XSdata('mat', groups)
    xsdata.order = 0
    xsdata.set_total([0.2])
    xsdata.set_absorption([0.1])
    xsdata.set_scatter_matrix(np.array([[[0.1]]]))
    mgxs_library = openmc.MGXSLibrary(groups)
    mgxs_library.add_xsdata(xsdata)
    mgxs_library.export_to_hdf5('mgxs.h5')
    mat = openmc.Material()
    mat.set_density('macro', 1.0)
    mat.add_macroscopic('mat')
    materials = openmc.Materials([mat])
    materials.cross_sections = 'mgxs.h5'

    width = 5.0
    inner = openmc.Sphere(r=3.098)
    outer = openmc.Sphere(r=3.1)
    box = openmc.model.RectangularParallelepiped(
        -width, width, -width, width, -width, width,
        boundary_type='reflective')
    shell = openmc.Cell(fill=mat, region=+inner & -outer)
    core = openmc.Cell(fill=mat, region=-inner)
    rest = openmc.Cell(fill=mat, region=+outer & -box)
    cells = [shell, core, rest]

    mesh = openmc.RegularMesh()
    mesh.dimension = (8, 8, 8)
    mesh.lower_left = (-width, -width, -width)
    mesh.upper_right = (width, width, width)

    settings = openmc.Settings()
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    settings.batches = 10
    settings.inactive = 5
    settings.particles = 50
    settings.random_ray['distance_active'] = 30.0
    settings.random_ray['distance_inactive'] = 30.0
    settings.random_ray['ray_source'] = openmc.IndependentSource(
        space=openmc.stats.Box(mesh.lower_left, mesh.upper_right))
    settings.random_ray['source_region_meshes'] = [(mesh, cells)]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(mesh.lower_left, mesh.upper_right),
        energy=openmc.stats.Discrete([1.0e6], [1.0]),
        constraints={'domains': [shell]},
        strength=1.0e18)
    settings.weight_window_generators = openmc.WeightWindowGenerator(
        method='fw_cadis', mesh=mesh, max_realizations=settings.batches)

    model = openmc.Model(openmc.Geometry(cells), materials, settings)
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
