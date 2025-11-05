import shutil
from pathlib import Path

import openmc
import openmc.lib
import pytest

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(), reason="DAGMC CAD geometry is not enabled."
)


@pytest.mark.parametrize("absolute", [True, False])
def test_model_h5m_in_subdirectory(run_in_tmpdir, request, absolute):
    # Create new subdirectory and copy h5m file there
    h5m = Path(request.fspath).parent / "dagmc.h5m"
    subdir = Path("h5m")
    subdir.mkdir()
    shutil.copy(h5m, subdir)

    # Create simple model with h5m file in subdirectory
    if absolute:
        dag_univ = openmc.DAGMCUniverse((subdir / "dagmc.h5m").absolute())
    else:
        dag_univ = openmc.DAGMCUniverse(subdir / "dagmc.h5m")
    model = openmc.Model()
    model.geometry = openmc.Geometry(dag_univ.bounded_universe())
    mat1 = openmc.Material(name="41")
    mat1.add_nuclide("H1", 1.0)
    mat2 = openmc.Material(name="no-void fuel")
    mat2.add_nuclide("U235", 1.0)
    model.materials = [mat1, mat2]
    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000

    # Make sure model can load
    model.export_to_model_xml()
    openmc.lib.init(["model.xml"])
    openmc.lib.finalize()
