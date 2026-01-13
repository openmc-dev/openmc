from pathlib import Path

import openmc


def test_statepoint_batches(run_in_tmpdir):
    # Create a minimal model
    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 4.5)
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 100

    # Specify when statepoints should be written
    model.settings.statepoint = {'batches': [3, 6, 9]}

    # Run model and ensure that statepoints are created
    model.run()
    sp_files = ['statepoint.03.h5', 'statepoint.06.h5', 'statepoint.09.h5']
    for f in sp_files:
        assert Path(f).is_file()
