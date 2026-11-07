import random
import xml.etree.ElementTree as ET

import numpy as np
import pytest
import openmc
import openmc.lib


SURFACE_IDS = [1, 2, 3, 4, 5, 6, 7]


def random_expression(rng, depth):
    """Random region expression using minimal parentheses so that operator
    precedence and complements of unparenthesized expressions are exercised"""
    if depth == 0 or rng.random() < 0.25:
        return f"{rng.choice(['', '-'])}{rng.choice(SURFACE_IDS)}"
    op = rng.choice([' ', ' | '])
    terms = [random_expression(rng, depth - 1) for _ in range(rng.randint(2, 3))]
    # A union only needs parentheses inside an intersection
    if op == ' ':
        terms = [f'({t})' if '|' in t else t for t in terms]
    expr = op.join(terms)
    if rng.random() < 0.3:
        expr = f'~({expr})'
    return expr


@pytest.mark.parametrize('seed', range(20))
def test_region_expression_parsing(run_in_tmpdir, seed):
    """Region expressions with minimal parentheses and complements are
    interpreted the same way by the C++ and Python parsers"""
    if seed == 0:
        # Complement of an expression mixing intersection and union
        expr = '~(1 3 | -5)'
    else:
        expr = random_expression(random.Random(seed), 4)

    surfaces = {s.id: s for s in [
        openmc.XPlane(-1.0, surface_id=1),
        openmc.XPlane(2.0, surface_id=2),
        openmc.YPlane(0.5, surface_id=3),
        openmc.YPlane(-2.5, surface_id=4),
        openmc.ZPlane(1.5, surface_id=5),
        openmc.Sphere(x0=1.0, r=3.0, surface_id=6),
        openmc.ZCylinder(y0=-1.0, r=2.0, surface_id=7),
        openmc.Sphere(r=10.0, boundary_type='vacuum', surface_id=99),
    ]}
    region = openmc.Region.from_expression(expr, surfaces)
    model = openmc.Model()
    model.geometry = openmc.Geometry([
        openmc.Cell(cell_id=1, region=-surfaces[99] & region),
        openmc.Cell(cell_id=2, region=-surfaces[99] & ~region),
    ])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1
    model.settings.batches = 1
    model.export_to_model_xml()

    # Give OpenMC the raw expressions rather than the fully parenthesized
    # strings written by the Python API
    tree = ET.parse('model.xml')
    for cell in tree.getroot().iter('cell'):
        if cell.get('id') == '1':
            cell.set('region', f'-99 {expr}')
        else:
            cell.set('region', f'-99 ~({expr})')
    tree.write('model.xml')

    points = np.random.default_rng(seed).uniform(-5.7, 5.7, size=(300, 3))
    openmc.lib.init()
    try:
        for p in points:
            cell, _ = openmc.lib.find_cell(p)
            expected = 1 if tuple(p) in region else 2
            assert cell.id == expected, f'{expr} at {p}'
    finally:
        openmc.lib.finalize()
