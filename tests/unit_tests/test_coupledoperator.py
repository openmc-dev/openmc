import openmc
import openmc.deplete


def test_diff_burnable_mats():
    """Tests the volumes assigned to the materials match the cell volumes"""
    mat1 = openmc.Material()
    mat1.add_element("Ag", 1, percent_type="ao")
    mat1.set_density("g/cm3", 10.49)
    mat1.depletable = True

    mat2 = openmc.Material()
    mat2.add_element("Ag", 1, percent_type="ao")
    mat2.set_density("g/cm3", 10.49)
    materials = openmc.Materials([mat1, mat2])

    sph1 = openmc.Sphere(r=1.0)
    sph2 = openmc.Sphere(r=2.0, x0=3)
    sph3 = openmc.Sphere(r=5.0, boundary_type="vacuum")

    cell1 = openmc.Cell(region=-sph1, cell_id=1)
    cell1.volume = 4.19
    cell2 = openmc.Cell(region=-sph2, cell_id=2)
    cell2.volume = 33.51
    cell3 = openmc.Cell(region=-sph3 & +sph1 & +sph2, cell_id=3)
    cell3.volume = 485.9

    cell1.fill = mat1
    cell2.fill = mat1
    cell3.fill = mat2

    geometry = openmc.Geometry([cell1, cell2, cell3])

    model = openmc.model.Model(geometry, materials)

    operator = openmc.deplete.CoupledOperator(
        model=model,
        normalization_mode="source-rate",
        reduce_chain=True,
        reduce_chain_level=5,
        diff_burnable_mats=True,
    )

    all_cells = operator.geometry.get_all_cells()
    assert all_cells[1].fill.volume == 4.19
    assert all_cells[2].fill.volume == 33.51
    # mat2 is not depletable
    assert all_cells[3].fill.volume == None
