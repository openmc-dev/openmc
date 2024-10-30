import pkg_resources
from pathlib import Path

import numpy as np
import pytest

import openmc
from openmc import ZPlane, YPlane, XPlane, Cell

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")

def test_model_differentiate_with_DAGMC():
    PITCH = 1.26

    mats = {}
    mats["no-void fuel"] = openmc.Material(1, "no-void fuel")
    mats["no-void fuel"].add_nuclide("U235", 0.03)
    mats["no-void fuel"].add_nuclide("U238", 0.97)
    mats["no-void fuel"].add_nuclide("O16", 2.0)
    mats["no-void fuel"].set_density("g/cm3", 10.0)
    mats["no-void fuel"].name = "Fuel"

    mats["41"] = openmc.Material(name="h2o")
    mats["41"].add_nuclide("H1", 2.0)
    mats["41"].add_element("O", 1.0)
    mats["41"].set_density("g/cm3", 1.0)
    mats["41"].add_s_alpha_beta("c_H_in_H2O")
    mats["41"].name = "41"

    p = pkg_resources.resource_filename(__name__, "dagmc.h5m")

    daguniv = openmc.DAGMCUniverse(p,auto_geom_ids=True,)

    def pattern(center, bc):
        bc_ = {
            "top": "transmission",
            "bottom": "transmission",
            "left": "transmission",
            "right": "transmission",
        }
        bc_ |= bc
        box = (
            -XPlane(center[0] + PITCH / 2, boundary_type=bc_["right"])
            & +XPlane(center[0] - PITCH / 2, boundary_type=bc_["left"])
            & -YPlane(center[1] + PITCH / 2, boundary_type=bc_["top"])
            & +YPlane(center[1] - PITCH / 2, boundary_type=bc_["bottom"])
            & -ZPlane(5, boundary_type="reflective")
            & +ZPlane(-5, boundary_type="reflective")
        )
        cell = Cell(region=box, fill=daguniv)
        cell.translation = [*center, 0]
        return [cell]

    root = openmc.Universe(
        cells=[
            *pattern((-PITCH / 2, -PITCH / 2),
                     bc={"left": "reflective", "bottom": "reflective"}),
            *pattern((-PITCH / 2, PITCH / 2),
                     bc={"left": "reflective", "top": "reflective"}),
            *pattern((PITCH / 2, PITCH / 2),
                     bc={"right": "reflective", "top": "reflective"}),
            *pattern((PITCH / 2, -PITCH / 2),
                     bc={"right": "reflective", "bottom": "reflective"}),
        ]
    )

    point = openmc.stats.Point((0, 0, 0))
    source = openmc.IndependentSource(space=point)

    settings = openmc.Settings()
    settings.source = source
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

    p = Path("differentiate_depletable_mats/divide_equally")
    p.mkdir(parents=True, exist_ok=True)
    model.init_lib()
    model.sync_dagmc_universe()
    model.calculate_volumes(cwd=p)
    volume_before = np.sum([m.volume for m in model.materials if m.name == "no-void fuel"])
    nmat = len(model.materials)
    model.differentiate_depletable_mats(diff_volume_method="divide equally")
    volume_after = np.sum([m.volume for m in model.materials if "fuel" in m.name])
    assert len(model.materials) == nmat + 3
    assert np.isclose(volume_before, volume_after)
    model.finalize_lib()
