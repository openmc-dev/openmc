"""Tally scoring in a void region while running in multi-group mode."""

import openmc
import pytest


# Speed of light in cm/s, used only as a physical upper bound on the neutron
# speed implied by a scored inverse velocity.
C_LIGHT = 2.99792458e10


@pytest.fixture
def one_group_lib():
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    xsdata = openmc.XSdata('slab_mat', groups)
    xsdata.order = 0
    xsdata.set_total([1.0])
    xsdata.set_absorption([0.5])
    xsdata.set_scatter_matrix([[[0.5]]])
    # Deliberately no set_inverse_velocity: the library then carries no
    # inverse-velocity data, so every region falls back to the approximate
    # group-average value in default_inverse_velocity_.

    mg_library = openmc.MGXSLibrary(groups)
    mg_library.add_xsdata(xsdata)
    name = 'mgxs.h5'
    mg_library.export_to_hdf5(name)
    yield name


@pytest.fixture
def void_model(one_group_lib):
    """Slab of material adjacent to a void region."""
    model = openmc.Model()
    mat = openmc.Material(name='slab_material')
    mat.set_density('macro', 1.0)
    mat.add_macroscopic('slab_mat')

    model.materials = openmc.Materials([mat])
    model.materials.cross_sections = one_group_lib

    x_min = openmc.XPlane(x0=-10.0, boundary_type='vacuum')
    x_mid = openmc.XPlane(x0=0.0)
    x_max = openmc.XPlane(x0=10.0, boundary_type='vacuum')
    y_min = openmc.YPlane(y0=-10.0, boundary_type='vacuum')
    y_max = openmc.YPlane(y0=10.0, boundary_type='vacuum')
    z_min = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')
    z_max = openmc.ZPlane(z0=10.0, boundary_type='vacuum')

    box = +y_min & -y_max & +z_min & -z_max
    solid = openmc.Cell(name='solid', fill=mat, region=+x_min & -x_mid & box)
    void = openmc.Cell(name='void', region=+x_mid & -x_max & box)
    model.geometry = openmc.Geometry([solid, void])

    model.settings = openmc.Settings()
    model.settings.energy_mode = 'multi-group'
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 5
    model.settings.particles = 100

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((-5.0, 0.0, 0.0))
    model.settings.source = source
    return model


def _cell(model, name):
    return next(c for c in model.geometry.get_all_cells().values()
                if c.name == name)


def test_flux_tally_in_void(run_in_tmpdir, void_model):
    """Flux is scored in a void cell instead of indexing absent material data."""
    tally = openmc.Tally()
    tally.filters = [openmc.CellFilter(_cell(void_model, 'void'))]
    tally.scores = ['flux']
    void_model.tallies = [tally]

    void_model.run(apply_tally_results=True)

    # Particles born in the slab stream into the void, so the flux there is
    # finite and positive rather than a crash or a zero-filled bin.
    assert tally.mean.squeeze() > 0.0


def test_inverse_velocity_in_void(run_in_tmpdir, void_model):
    """Inverse velocity is defined in a void and scores there.

    It depends only on the particle's speed, not on any material, so a void
    region is not a reason to drop it. Before this was handled, scoring it
    indexed the macroscopic cross section table with MATERIAL_VOID and the run
    aborted.
    """
    void_tally = openmc.Tally(name='void')
    void_tally.filters = [openmc.CellFilter(_cell(void_model, 'void'))]
    void_tally.scores = ['flux', 'inverse-velocity']

    solid_tally = openmc.Tally(name='solid')
    solid_tally.filters = [openmc.CellFilter(_cell(void_model, 'solid'))]
    solid_tally.scores = ['flux', 'inverse-velocity']

    void_model.tallies = [void_tally, solid_tally]
    void_model.run(apply_tally_results=True)

    void_flux, void_inv_v = void_tally.mean.ravel()
    solid_flux, solid_inv_v = solid_tally.mean.ravel()

    # The score is present, not silently zeroed.
    assert void_inv_v > 0.0

    # Each track contributes its length times the same group constant, so the
    # inverse velocity per unit flux is that constant exactly. The library
    # carries no inverse-velocity data, so the material falls back to the same
    # default the void uses and the two must agree to round-off.
    assert (void_inv_v / void_flux ==
            pytest.approx(solid_inv_v / solid_flux, rel=1e-9))

    # That constant is a real inverse velocity: the speed it implies is
    # positive and below the speed of light.
    speed = void_flux / void_inv_v
    assert 0.0 < speed < C_LIGHT


def test_events_scored_in_void(run_in_tmpdir, void_model):
    """Counting scoring events needs no material, so a void cell still counts.

    'events' is score = 1.0 per scoring event in both energy modes; nothing
    about it refers to the material.
    """
    tally = openmc.Tally()
    tally.filters = [openmc.CellFilter(_cell(void_model, 'void'))]
    tally.scores = ['events']
    void_model.tallies = [tally]

    void_model.run(apply_tally_results=True)

    assert tally.mean.squeeze() > 0.0


def test_material_dependent_scores_in_void(run_in_tmpdir, void_model):
    """Scores that need material data are zero in a void cell."""
    tally = openmc.Tally()
    tally.filters = [openmc.CellFilter(_cell(void_model, 'void'))]
    tally.scores = ['total', 'absorption']
    void_model.tallies = [tally]

    void_model.run(apply_tally_results=True)

    assert (tally.mean.squeeze() == 0.0).all()
