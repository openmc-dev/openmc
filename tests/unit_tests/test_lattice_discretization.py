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

    assert rlat_clone.get_universe((0, 2)) == rlat_clone.get_universe((2, 2))
    assert rlat_clone.get_universe((0, 2)) != rlat_clone.get_universe((0, 0))
    assert rlat_clone.get_universe((0, 2)) != rlat_clone.get_universe((0, 1))
    assert rlat_clone.get_universe((0, 2)) != rlat_clone.get_universe((1, 0))
    assert rlat_clone.get_universe((0, 2)) != rlat_clone.get_universe((1, 1))
    assert rlat_clone.get_universe((0, 2)) != rlat_clone.get_universe((1, 2))
    assert rlat_clone.get_universe((0, 2)) != rlat_clone.get_universe((2, 0))
    assert rlat_clone.get_universe((0, 2)) != rlat_clone.get_universe((2, 1))

def test_discretization_lns_using_names(rlat2):

    rlat_clone = rlat2.clone()
    rlat_clone.discretize(strategy="lns", key=lambda univ: univ.name)

    assert rlat_clone.get_universe((0, 2)) == rlat_clone.get_universe((2, 2))

    rlat_clone = rlat2.clone()
    rlat_clone.get_universe((0, 1)).name="u1"
    rlat_clone.get_universe((1, 0)).name="u1"
    rlat_clone.discretize(strategy="lns", key=lambda univ: univ.name)

    assert rlat_clone.get_universe((0, 2)) == rlat_clone.get_universe((2, 0))
    assert rlat_clone.get_universe((2, 2)) == rlat_clone.get_universe((0, 0))
    assert rlat_clone.get_universe((0, 2)) == rlat_clone.get_universe((2, 2))
    assert rlat_clone.get_universe((1, 2)) == rlat_clone.get_universe((1, 0))
    assert rlat_clone.get_universe((2, 1)) == rlat_clone.get_universe((0, 1))
    assert rlat_clone.get_universe((1, 2)) == rlat_clone.get_universe((0, 1))
    assert rlat_clone.get_universe((0, 0)) != rlat_clone.get_universe((0, 1))

def test_discretization_lns_with_neighbor_list(rlat2):

    rlat_clone = rlat2.clone()
    u1 = rlat_clone.get_universe((1, 0))
    u2 = rlat_clone.get_universe((0, 0))
    rlat_clone.discretize(strategy="lns",
         lattice_neighbors=[u1,u1,u2,u1,u2,u2,u1,u2])

    assert rlat_clone.get_universe((0, 1)) == rlat_clone.get_universe((0, 0))
    assert rlat_clone.get_universe((1, 1)) != rlat_clone.get_universe((0, 0))

