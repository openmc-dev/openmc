import openmc
import openmc.stats
import h5py


def test_musurface(run_in_tmpdir):
    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=None)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 100
    model.settings.batches = 1
    E = 1.0
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.delta_function(E),
        weight=openmc.stats.delta_function(100)
    )
    model.settings.run_mode = "fixed source"
    model.settings.surf_source_write = {
        "max_particles": 100,
    }

    # Run OpenMC
    sp_filename = model.run()

    # All contributions should show up in last bin
    with h5py.File("surface_source.h5", "r") as f:
        source = f["source_bank"]

        assert len(source) == 100

        for point in source:
            assert point["wgt"] == 100.0


