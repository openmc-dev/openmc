from math import sqrt
import xml.etree.ElementTree as ET

import openmc
import pytest

from tests.unit_tests.test_lattice import zr, pincell1, pincell2, rlat2


def test_discretization_degenerate(rlat2):

    rlat_clone = rlat2.clone()
    rlat_clone.discretize()

    assert rlat_clone.get_universe((0, 0)).id !=\
         rlat_clone.get_universe((1, 0)).id
    assert rlat_clone.get_universe((1, 0)).id !=\
         rlat_clone.get_universe((0, 1)).id
    assert rlat_clone.get_universe((0, 1)).id !=\
         rlat_clone.get_universe((0, 0)).id
    assert rlat_clone.get_universe((0, 2)).id !=\
         rlat_clone.get_universe((2, 0)).id
    assert rlat_clone.get_universe((2, 1)).id !=\
         rlat_clone.get_universe((1, 2)).id


def test_discretization_skip_universe(rlat2):

    rlat_clone = rlat2.clone()
    rlat_clone.discretize(
         universes_to_ignore=[rlat_clone.get_universe((0, 0))])

    assert rlat_clone.get_universe((0, 0)) == rlat_clone.get_universe((0, 1))
    assert rlat_clone.get_universe((0, 1)) == rlat_clone.get_universe((2, 1))
    assert rlat_clone.get_universe((0, 1)) == rlat_clone.get_universe((1, 2))


def test_discretization_clone_only_some_materials(rlat2):

    rlat_clone = rlat2.clone()
    fuel1 = next(iter(rlat_clone.get_universe((0, 0)).cells.values())).fill
    rlat_clone.discretize(materials_to_clone=[fuel1])

    assert next(reversed(rlat_clone.get_universe((0, 0)).cells.values())).fill\
         == next(reversed(rlat_clone.get_universe((1, 0)).cells.values())).fill
    assert next(iter(rlat_clone.get_universe((0, 0)).cells.values())).fill\
         != next(iter(rlat_clone.get_universe((1, 0)).cells.values())).fill


def test_discretization_lns(rlat2):

    rlat_clone = rlat2.clone()
    rlat_clone.discretize(strategy="lns")

    assert rlat_clone.get_universe((0, 2)) == \
         next(iter(rlat_clone.get_universe((2, 2)).cells.values())).fill


def test_discretization_lns_using_names(rlat2):

    rlat_clone = rlat2.clone()
    rlat_clone.discretize(strategy="lns", attribute="name")

    assert rlat_clone.get_universe((0, 2)) == \
         rlat_clone.get_universe((2, 2))

    rlat_clone = rlat2.clone()
    rlat_clone.get_universe((1, 0)).name="u1"
    rlat_clone.discretize(strategy="lns", attribute="name")

    assert rlat_clone.get_universe((0, 2)) == \
         next(iter(rlat_clone.get_universe((2, 2)).cells.values())).fill


def test_discretization_lns_with_neighbor_list(rlat2):

    rlat_clone = rlat2.clone()
    u1 = rlat_clone.get_universe((1, 0))
    u2 = rlat_clone.get_universe((0, 0))
    rlat_clone.discretize(strategy="lns",
         lattice_neighbors=[u1.id,u1.id,u2.id,u1.id,u2.id,u2.id,u1.id,u2.id])

    assert rlat_clone.get_universe((0, 1)) == \
         next(iter(rlat_clone.get_universe((0, 0)).cells.values())).fill
    assert rlat_clone.get_universe((1, 1)) != \
         next(iter(rlat_clone.get_universe((0, 0)).cells.values())).fill


def test_discretization_lns_nested_lattices(rlat2):

    # Name a universe in rlat2 to be able to compare
    rlat2.get_universe((1, 2)).name = "universe"

    # Create a lattice of lattices
    main_lattice = openmc.RectLattice(name="main lattice")
    u1 = openmc.Universe()
    cell = openmc.Cell()
    cell.fill=rlat2
    u1.add_cell(cell)
    u2 = u1.clone()
    main_lattice.pitch = (20, 20)
    main_lattice.universes = [[u1,u1,u1], [u2,u1,u2], [u1,u1,u2]]

    main_lattice.discretize(strategy="lns")

    # Both involved a transposition from (2,0) neighbor pattern
    assert next(iter(main_lattice.get_universe((0, 0)).cells.values())).fill\
         == next(iter(main_lattice.get_universe((2, 2)).cells.values())).fill

    # Transposition creates a clone so the two universes cant match
    assert main_lattice.get_universe((2, 0)) != \
         next(iter(main_lattice.get_universe((2, 2)).cells.values())).fill

    # Different neighbor pattern
    assert main_lattice.get_universe((1, 1)) != \
         next(iter(main_lattice.get_universe((0, 0)).cells.values())).fill

    # Look inside inner lattices
    for i in range(3):
        for j in range(3):
            assert next(iter(main_lattice.get_universe((0, 2)).cells.values())\
                 ).fill.get_universe((i, j)).name == \
                 next(iter(next(iter(main_lattice.get_universe((2, 2)).cells\
                 .values())).fill.cells.values())).fill.get_universe((2-j,
                 2-i)).name

