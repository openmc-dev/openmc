import pkg_resources

import numpy as np
import pytest

import openmc

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")

@pytest.fixture()
def model():
    PITCH = 1.26

    mats = {}
    mats["no-void fuel"] = openmc.Material(1, name="no-void fuel")
    mats["no-void fuel"].add_nuclide("U235", 0.03)
    mats["no-void fuel"].add_nuclide("U238", 0.97)
    mats["no-void fuel"].add_nuclide("O16", 2.0)
    mats["no-void fuel"].set_density("g/cm3", 10.0)

    mats["41"] = openmc.Material(name="41")
    mats["41"].add_nuclide("H1", 2.0)
    mats["41"].add_element("O", 1.0)
    mats["41"].set_density("g/cm3", 1.0)
    mats["41"].add_s_alpha_beta("c_H_in_H2O")

    p = pkg_resources.resource_filename(__name__, "dagmc.h5m")

    daguniv = openmc.DAGMCUniverse(p,auto_geom_ids=True,)

    lattice = openmc.RectLattice()
    lattice.dimension = [2, 2]
    lattice.lower_left = [-PITCH, -PITCH]
    lattice.pitch = [PITCH, PITCH]
    lattice.universes = [
        [daguniv, daguniv],
        [daguniv, daguniv]]

    box = openmc.model.RectangularParallelepiped(-PITCH, PITCH, -PITCH, PITCH, -5, 5)

    root = openmc.Universe(cells=[openmc.Cell(region= -box, fill=lattice)])

    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 10
    settings.particles = 1000

    ll, ur = root.bounding_box
    mat_vol = openmc.VolumeCalculation([mats["no-void fuel"]], 1000000, ll, ur)
    cell_vol = openmc.VolumeCalculation(
        list(root.cells.values()), 1000000, ll, ur)
    settings.volume_calculations = [mat_vol, cell_vol]

    model = openmc.Model()
    model.materials = openmc.Materials(mats.values())
    model.geometry = openmc.Geometry(root=root)
    model.settings = settings
    return model


def test_model_differentiate_depletable_with_dagmc(model, run_in_tmpdir):
    try:
        model.init_lib()
        model.sync_dagmc_universes()
        model.calculate_volumes()

        # Get the volume of the no-void fuel material before differentiation
        volume_before = np.sum([m.volume for m in model.materials if m.name == "no-void fuel"])

        # Differentiate the depletable materials
        model.differentiate_depletable_mats(diff_volume_method="divide equally")
        # Get the volume of the no-void fuel material after differentiation
        volume_after = np.sum([m.volume for m in model.materials if "fuel" in m.name])
        assert np.isclose(volume_before, volume_after)
        assert len(model.materials) == 4*2 +1
    finally:
        model.finalize_lib()


def test_model_differentiate_with_dagmc(model, run_in_tmpdir):
    root = model.geometry.root_universe
    ll, ur = root.bounding_box
    try:
        model.init_lib()
        model.sync_dagmc_universes()
        model.calculate_volumes()
        # Get the volume of the no-void fuel material before differentiation
        volume_before = np.sum([m.volume for m in model.materials if m.name == "no-void fuel"])

        # Differentiate all the materials
        model.differentiate_mats(depletable_only=False)

        # Get the volume of the no-void fuel material after differentiation
        mat_list = [m for m in model.materials]
        mat_vol = openmc.VolumeCalculation(mat_list, 1000000, ll, ur)
        cell_vol = openmc.VolumeCalculation(list(root.cells.values()), 1000000, ll, ur)
        model.settings.volume_calculations = [mat_vol, cell_vol]
        model.init_lib() # need to reinitialize the lib after differentiating the materials
        model.calculate_volumes()
        volume_after = np.sum([m.volume for m in model.materials if "fuel" in m.name])
        assert np.isclose(volume_before, volume_after)
        assert len(model.materials) == 4*2 + 4
    finally:
        model.finalize_lib()


def test_bad_override_cell_id(model, run_in_tmpdir):
    try:
        model.init_lib()
        model.sync_dagmc_universes()
        for univ in model.geometry.get_all_universes().values():
            if isinstance(univ, openmc.DAGMCUniverse):
                break
        with pytest.raises(ValueError, match="Cell ID '1' not found in DAGMC universe"):
            univ.material_overrides = {1 : model.materials[0]}
    finally:
        model.finalize_lib()