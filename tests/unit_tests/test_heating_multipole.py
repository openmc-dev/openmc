import openmc
import pytest


@pytest.mark.parametrize("multipole", [True, False])
def test_heating_multipole_warning(run_in_tmpdir, capfd, mpi_intracomm,
                                   multipole):
    # Heating is not scored at energies where cross sections are evaluated
    # from windowed multipole data, so heating scores should give a warning
    # when multipole data is in use
    mat = openmc.Material()
    mat.add_nuclide("U235", 1.0)
    mat.set_density("g/cm3", 10.0)
    sphere = openmc.Sphere(r=10.0, boundary_type="vacuum")
    model = openmc.Model()
    model.geometry = openmc.Geometry([openmc.Cell(fill=mat, region=-sphere)])
    model.settings.particles = 100
    model.settings.batches = 1
    model.settings.temperature = {"multipole": multipole}

    tally = openmc.Tally()
    tally.scores = ["heating", "heating-local", "fission"]
    model.tallies = [tally]

    model.init_lib(output=False, intracomm=mpi_intracomm)
    model.finalize_lib()

    err = capfd.readouterr().err
    assert ("'heating' score" in err) == multipole
    assert ("'heating-local' score" in err) == multipole
    assert "'fission' score" not in err
