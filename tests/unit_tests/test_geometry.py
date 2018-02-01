import xml.etree.ElementTree as ET

import openmc
import pytest


def test_volume(uo2):
    """Test adding volume information from a volume calculation."""
    # Create model with nested spheres
    model = openmc.model.Model()
    model.materials.append(uo2)
    inner = openmc.Sphere(R=1.)
    outer = openmc.Sphere(R=2., boundary_type='vacuum')
    c1 = openmc.Cell(fill=uo2, region=-inner)
    c2 = openmc.Cell(region=+inner & -outer)
    u = openmc.Universe(cells=[c1, c2])
    model.geometry.root_universe = u
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.Source(space=openmc.stats.Point())

    ll, ur = model.geometry.bounding_box
    assert ll == pytest.approx((-outer.r, -outer.r, -outer.r))
    assert ur == pytest.approx((outer.r, outer.r, outer.r))
    model.settings.volume_calculations

    for domain in (c1, uo2, u):
        # Run stochastic volume calculation
        volume_calc = openmc.VolumeCalculation(
            domains=[domain], samples=1000, lower_left=ll, upper_right=ur)
        model.settings.volume_calculations = [volume_calc]
        model.export_to_xml()
        openmc.calculate_volumes()

        # Load results and add volume information
        volume_calc.load_results('volume_1.h5')
        model.geometry.add_volume_information(volume_calc)

        # get_nuclide_densities relies on volume information
        nucs = set(domain.get_nuclide_densities())
        assert not (nucs ^ {'U235', 'O16'})


def test_export_xml():
    pass


def test_find():
    pass


def test_get_instances():
    pass


def test_get_all_universes():
    pass


def test_get_all_materials():
    pass


def test_get_all_material_cells():
    pass


def test_get_all_material_universes():
    pass


def test_get_all_lattices():
    pass


def test_get_all_surfaces():
    pass


def test_get_materials_by_name():
    pass


def test_get_cells_by_name():
    pass


def test_get_cells_by_fill_name():
    pass


def test_get_universes_by_name():
    pass


def test_get_lattices_by_name():
    pass


def test_clone():
    pass


def test_determine_paths(cell_with_lattice):
    cells, mats, univ, lattice = cell_with_lattice
    u = openmc.Universe(cells=[cells[-1]])
    geom = openmc.Geometry(u)

    geom.determine_paths()
    assert len(cells[0].paths) == 4
    assert len(cells[1].paths) == 4
    assert len(cells[2].paths) == 1
    assert len(mats[0].paths) == 1
    assert len(mats[-1].paths) == 4
