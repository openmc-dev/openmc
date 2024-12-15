import lxml.etree as ET
import numpy as np
import openmc
import pytest

from tests.unit_tests import assert_unbounded


def test_basic():
    c1 = openmc.Cell()
    c2 = openmc.Cell()
    c3 = openmc.Cell()
    u = openmc.Universe(name="cool", cells=(c1, c2, c3))
    assert u.name == "cool"

    cells = set(u.cells.values())
    assert not (cells ^ {c1, c2, c3})

    # Test __repr__
    repr(u)

    with pytest.raises(TypeError):
        u.add_cell(openmc.Material())
    with pytest.raises(TypeError):
        u.add_cells(c1)

    u.remove_cell(c3)
    cells = set(u.cells.values())
    assert not (cells ^ {c1, c2})

    u.clear_cells()
    assert not set(u.cells)


def test_bounding_box():
    cyl1 = openmc.ZCylinder(r=1.0)
    cyl2 = openmc.ZCylinder(r=2.0)
    c1 = openmc.Cell(region=-cyl1)
    c2 = openmc.Cell(region=+cyl1 & -cyl2)

    u = openmc.Universe(cells=[c1, c2])
    ll, ur = u.bounding_box
    assert ll == pytest.approx((-2.0, -2.0, -np.inf))
    assert ur == pytest.approx((2.0, 2.0, np.inf))

    u = openmc.Universe()
    assert_unbounded(u)


def test_plot(run_in_tmpdir, sphere_model):

    # model with -inf and inf in the bounding box
    pincell = openmc.examples.pwr_pin_cell()
    materials = pincell.materials

    mat_colors = {
        materials[0]: (200, 1, 1),
        materials[1]: "gray",
        materials[2]: "limegreen",
    }

    for basis in ("xy", "yz", "xz"):
        plot = pincell.geometry.root_universe.plot(
            colors=mat_colors,
            color_by="material",
            legend=True,
            pixels=(10, 10),
            basis=basis,
            outline=True,
            axis_units="m",
        )
        assert plot.xaxis.get_label().get_text() == f"{basis[0]} [m]"
        assert plot.yaxis.get_label().get_text() == f"{basis[1]} [m]"

    # model with no inf values in bounding box
    m = sphere_model.materials[0]
    univ = sphere_model.geometry.root_universe

    colors = {m: "limegreen"}

    for basis in ("xy", "yz", "xz"):
        plot = univ.plot(
            colors=colors,
            color_by="cell",
            legend=False,
            pixels=100,
            basis=basis,
            outline=False,
        )
        assert plot.xaxis.get_label().get_text() == f"{basis[0]} [cm]"
        assert plot.yaxis.get_label().get_text() == f"{basis[1]} [cm]"

    msg = "Must pass 'colors' dictionary if you are adding a legend via legend=True."
    # This plot call should fail as legend is True but colors is None
    with pytest.raises(ValueError, match=msg):
        univ.plot(
            color_by="cell",
            legend=True,
            pixels=100,
        )


def test_get_nuclides(uo2):
    c = openmc.Cell(fill=uo2)
    univ = openmc.Universe(cells=[c])
    nucs = univ.get_nuclides()
    assert nucs == ["U235", "O16"]


def test_cells():
    cells = [openmc.Cell() for i in range(5)]
    cells2 = [openmc.Cell() for i in range(3)]
    cells[0].fill = openmc.Universe(cells=cells2)
    u = openmc.Universe(cells=cells)
    assert not (set(u.cells.values()) ^ set(cells))

    all_cells = set(u.get_all_cells().values())
    assert not (all_cells ^ set(cells + cells2))


def test_get_all_materials(cell_with_lattice):
    cells, mats, univ, lattice = cell_with_lattice
    test_mats = set(univ.get_all_materials().values())
    assert not (test_mats ^ set(mats))


def test_get_all_universes():
    c1 = openmc.Cell()
    u1 = openmc.Universe(cells=[c1])
    c2 = openmc.Cell()
    u2 = openmc.Universe(cells=[c2])
    c3 = openmc.Cell(fill=u1)
    c4 = openmc.Cell(fill=u2)
    u3 = openmc.DAGMCUniverse(filename="")
    c5 = openmc.Cell(fill=u3)
    u4 = openmc.Universe(cells=[c3, c4, c5])

    univs = set(u4.get_all_universes().values())
    assert not (univs ^ {u1, u2, u3})


def test_clone():

    c1 = openmc.Cell(cell_id=1)
    c1.region = -openmc.ZCylinder(r=1.0)
    c2 = openmc.Cell(cell_id=2)
    c2.fill = openmc.Material()
    c3 = openmc.Cell()
    u1 = openmc.Universe(name="cool", cells=(c1, c2, c3))
    u1.volume = 1.0

    u2 = u1.clone()
    assert u2.name == u1.name
    assert u2.cells != u1.cells
    assert u2.get_all_materials() != u1.get_all_materials()
    assert u2.volume == u1.volume

    u2 = u1.clone(clone_materials=False)
    assert u2.get_all_materials() == u1.get_all_materials()

    u3 = u1.clone(clone_regions=False)
    assert next(iter(u3.cells.values())).region == next(iter(u1.cells.values())).region

    # Change attributes, make sure clone stays intact
    u1.volume = 2.0
    u1.name = "different name"
    assert u3.volume != u1.volume
    assert u3.name != u1.name

    # Test cloning a DAGMC universe
    dagmc_u = openmc.DAGMCUniverse(filename="", name="DAGMC universe")
    dagmc_u.volume = 1.0
    dagmc_u.auto_geom_ids = True
    dagmc_u.auto_mat_ids = True
    dagmc_u1 = dagmc_u.clone()
    assert dagmc_u1.name == dagmc_u.name
    assert dagmc_u1.volume == dagmc_u.volume
    assert dagmc_u1.auto_geom_ids == dagmc_u.auto_geom_ids
    assert dagmc_u1.auto_mat_ids == dagmc_u.auto_mat_ids

    # Change attributes, check the clone remained intact
    dagmc_u.name = "another name"
    dagmc_u.auto_geom_ids = False
    dagmc_u.auto_mat_ids = False
    dagmc_u.volume = 2.0
    assert dagmc_u1.name != dagmc_u.name
    assert dagmc_u1.volume != dagmc_u.volume
    assert dagmc_u1.auto_geom_ids != dagmc_u.auto_geom_ids
    assert dagmc_u1.auto_mat_ids != dagmc_u.auto_mat_ids


def test_create_xml(cell_with_lattice):
    cells = [openmc.Cell() for i in range(5)]
    u = openmc.Universe(cells=cells)

    geom = ET.Element("geom")
    u.create_xml_subelement(geom)
    cell_elems = geom.findall("cell")
    assert len(cell_elems) == len(cells)
    assert all(c.get("universe") == str(u.id) for c in cell_elems)
    assert not (set(c.get("id") for c in cell_elems) ^ set(str(c.id) for c in cells))


def test_get_nuclide_densities():
    surf = openmc.Sphere()
    material = openmc.Material()
    material.add_elements_from_formula("H2O")
    material.set_density("g/cm3", 1)
    cell = openmc.Cell(region=-surf, fill=material)
    universe = openmc.Universe(cells=[cell])
    with pytest.raises(RuntimeError):
        universe.get_nuclide_densities()
