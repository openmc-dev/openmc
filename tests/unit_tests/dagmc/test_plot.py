import pytest
import openmc


pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(), reason="DAGMC CAD geometry is not enabled."
)

def test_ploting_dagmc_geometry():
    """Test plotting a DAGMC geometry with OpenMC. This is different to CSG
    geometry plotting as the path to the DAGMC file needs handling."""

    with openmc.config.patch('resolve_paths', True):
        dag_universe = openmc.DAGMCUniverse(filename='dagmc.h5m')
    csg_with_dag_inside = dag_universe.bounded_universe()
    my_geometry = openmc.Geometry(csg_with_dag_inside)

    all_materials = []
    for mat_name in dag_universe.material_names:
        material = openmc.Material(name=mat_name)
        material.add_nuclide("Fe56", 1, "ao")
        material.set_density("g/cm3", 7)
        all_materials.append(material)
    my_materials = openmc.Materials(all_materials)

    my_source = openmc.IndependentSource()
    # putting the source at the center of the bounding box of the DAGMC
    my_source.space = openmc.stats.Point(dag_universe.bounding_box.center)

    my_settings = openmc.Settings(
        batches=10,
        particles=50,
        source=my_source
    )

    my_model = openmc.Model(
        geometry=my_geometry,
        materials=my_materials,
        settings=my_settings
    )

    my_model.plot()
