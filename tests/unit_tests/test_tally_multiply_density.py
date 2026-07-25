import numpy as np
import openmc
import pytest


def test_micro_macro_compare(run_in_tmpdir):
    # Create simple sphere model with H1 and H2
    mat = openmc.Material()
    mat.add_components({'H1': 1.0, 'H2': 1.0})
    mat.set_density('g/cm3', 1.0)
    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1000
    model.settings.batches = 10

    # Set up two reaction rate tallies, one that multplies by density and the
    # other that doesn't
    tally_macro = openmc.Tally()
    tally_macro.nuclides = ['H1', 'H2', 'H3']
    tally_macro.scores = ['total', 'elastic']
    tally_micro = openmc.Tally()
    tally_micro.nuclides = ['H1', 'H2', 'H3']
    tally_micro.scores = ['total', 'elastic']
    tally_micro.multiply_density = False
    model.tallies = [tally_macro, tally_micro]

    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        tally_macro = sp.tallies[tally_macro.id]
        tally_micro = sp.tallies[tally_micro.id]

    # Make sure multply_density attribute from statepoint is set correctly
    assert tally_macro.multiply_density
    assert not tally_micro.multiply_density

    # Dividing macro by density should give micro
    density = mat.get_nuclide_atom_densities()
    for nuc in ('H1', 'H2'):
        micro_derived = tally_macro.get_values(nuclides=[nuc]) / density[nuc]
        micro = tally_micro.get_values(nuclides=[nuc])
        assert micro_derived == pytest.approx(micro)

    # For macro tally, H3 scores should be zero
    assert np.all(tally_macro.get_values(nuclides=['H3']) == 0.0)

    # For micro tally, H3 scores should be positive
    assert np.all(tally_micro.get_values(nuclides=['H3']) > 0.0)


def test_overlay_material(run_in_tmpdir):
    """An overlay material bin should give the material's macroscopic rate."""
    # The geometry is filled with iron and contains no hydrogen at all, so the
    # overlay material can only be scored if it is decoupled from the geometry
    iron = openmc.Material()
    iron.add_nuclide('Fe56', 1.0)
    iron.set_density('g/cm3', 7.874)

    water = openmc.Material()
    water.add_components({'H1': 2.0, 'O16': 1.0})
    water.set_density('g/cm3', 1.0)

    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=iron, region=-sph)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1000
    model.settings.batches = 10

    # One tally with a single overlay material bin and one with a bin per
    # nuclide of that material. Summing the latter weighted by atom density
    # must reproduce the former.
    tally_mat = openmc.Tally()
    tally_mat.nuclides = [water]
    tally_mat.scores = ['total', 'absorption']
    tally_mat.multiply_density = False

    tally_nucs = openmc.Tally()
    tally_nucs.nuclides = ['H1', 'O16']
    tally_nucs.scores = ['total', 'absorption']
    tally_nucs.multiply_density = False

    model.tallies = [tally_mat, tally_nucs]

    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        tally_mat = sp.tallies[tally_mat.id]
        tally_nucs = sp.tallies[tally_nucs.id]

    # The overlay material is not in the geometry but must still be written
    assert water.id in {m.id for m in model._materials_for_export()}

    # The nuclide axis label round trips through the statepoint
    assert tally_mat.nuclides == [f'material:{water.id}']

    # The tolerance is loose because the atom densities used here come from
    # openmc.data.atomic_mass whereas the solver derives them from the AWR
    # values in the nuclear data files, which differ by a few parts in 1e5
    density = water.get_nuclide_atom_densities()
    expected = sum(density[nuc] * tally_nucs.get_values(nuclides=[nuc])
                   for nuc in ('H1', 'O16'))
    assert tally_mat.get_values() == pytest.approx(expected, rel=1e-4)

    # The scores are nonzero, i.e. the comparison above is not 0 == 0
    assert np.all(tally_mat.get_values() > 0.0)


def test_overlay_material_errors(run_in_tmpdir):
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    mat.set_density('g/cm3', 1.0)

    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 100
    model.settings.batches = 2

    # An overlay material requires multiply_density to be False
    tally = openmc.Tally()
    tally.nuclides = [mat]
    tally.scores = ['absorption']
    model.tallies = [tally]
    with pytest.raises(ValueError, match='multiply_density'):
        model.export_to_model_xml()

    # 'flux' does not depend on the nuclide
    tally.multiply_density = False
    tally.scores = ['flux']
    with pytest.raises(RuntimeError, match='flux'):
        model.run()

    # The overlay is only scored on the tracklength and collision paths
    tally.scores = ['absorption']
    tally.estimator = 'analog'
    with pytest.raises(RuntimeError, match='analog estimator'):
        model.run()
