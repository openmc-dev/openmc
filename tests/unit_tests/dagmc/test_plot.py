import pytest
import openmc
import openmc.lib


pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(), reason="DAGMC CAD geometry is not enabled."
)

def test_plotting_dagmc_geometry(request):
    """Test plotting a DAGMC geometry with OpenMC. This is different to CSG
    geometry plotting as the path to the DAGMC file needs handling."""

    dag_universe = openmc.DAGMCUniverse(request.path.parent / 'dagmc.h5m')
    csg_with_dag_inside = dag_universe.bounded_universe()
    model = openmc.Model()
    model.geometry = openmc.Geometry(csg_with_dag_inside)

    for mat_name in dag_universe.material_names:
        material = openmc.Material(name=mat_name)
        material.add_nuclide("Fe56", 1.0)
        material.set_density("g/cm3", 7.0)
        model.materials.append(material)

    # putting the source at the center of the bounding box of the DAGMC
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(dag_universe.bounding_box.center)
    )
    model.settings.batches = 10
    model.settings.particles = 50

    model.plot()
