import pytest
import openmc
import openmc.lib


pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(), reason="DAGMC CAD geometry is not enabled."
)

def test_plotting_dagmc_model(request):
    """Test plotting a DAGMC model with OpenMC. This is different to CSG
    model plotting as the path to the DAGMC file needs handling."""

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


def test_plotting_dagmc_universe(request):
    """Test plotting a DAGMCUniverse with OpenMC. This is different to plotting
    UniverseBase as the materials are not defined withing the DAGMCUniverse."""

    dag_universe = openmc.DAGMCUniverse(request.path.parent / 'dagmc.h5m')
    dag_universe.plot()


def test_plotting_geometry_filled_with_dagmc_universe(request):
    """Test plotting a geometry with OpenMC. This is an edge case when plotting
    geometry as often geometry objects don't include a DAGMCUniverse. The
    inclusion of a DAGMCUniverse requires special handling for the materials."""

    dag_universe = openmc.DAGMCUniverse(request.path.parent / 'dagmc.h5m', auto_geom_ids=True)

    sphere1 = openmc.Sphere(r=50.0)
    sphere2 = openmc.Sphere(r=60.0, boundary_type='vacuum')

    csg_material = openmc.Material(name='csg_material')
    csg_material.add_nuclide("H1", 1.0)

    cell1 = openmc.Cell(fill=dag_universe, region=-sphere1)
    cell2 = openmc.Cell(fill=dag_universe, region=+sphere1 & -sphere2)

    geometry = openmc.Geometry([cell1, cell2])

    geometry.plot()
