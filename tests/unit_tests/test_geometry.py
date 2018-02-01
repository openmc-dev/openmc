import xml.etree.ElementTree as ET

import openmc
import pytest


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
