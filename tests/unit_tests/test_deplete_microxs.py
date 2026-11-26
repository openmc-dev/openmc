"""Basic unit tests for openmc.deplete.IndependentOperator instantiation

Modifies and resets environment variable OPENMC_CROSS_SECTIONS
to a custom file with new depletion_chain node
"""

from os import remove
from pathlib import Path
from unittest.mock import patch

from lxml import etree as ET
import pytest
import openmc
from openmc.deplete import MicroXS, get_microxs_and_flux
from openmc.mixin import IDManagerMixin, IDWarning
import numpy as np

ONE_GROUP_XS = Path(__file__).parents[1] / "micro_xs_simple.csv"
CHAIN_FILE = Path(__file__).parents[1] / "chain_simple.xml"


def test_from_array():
    nuclides = [
        'U234',
        'U235',
        'U238',
        'U236',
        'O16',
        'O17',
        'I135',
        'Xe135',
        'Xe136',
        'Cs135',
        'Gd157',
        'Gd156']
    reactions = ['fission', '(n,gamma)']
    # These values are placeholders and are not at all
    # physically meaningful.
    data = np.array([[0.1, 0.],
                     [0.1, 0.],
                     [0.9, 0.],
                     [0.4, 0.],
                     [0., 0.],
                     [0., 0.],
                     [0., 0.1],
                     [0., 0.9],
                     [0., 0.],
                     [0., 0.],
                     [0., 0.1],
                     [0., 0.1]])
    data = data.reshape(12, 2, 1)

    MicroXS(data, nuclides, reactions)
    with pytest.raises(ValueError, match='Data array must be 3D'):
        MicroXS(data[:, 0], nuclides, reactions)


def test_csv():
    ref_xs = MicroXS.from_csv(ONE_GROUP_XS)
    ref_xs.to_csv('temp_xs.csv')
    temp_xs = MicroXS.from_csv('temp_xs.csv')
    assert np.all(ref_xs.data == temp_xs.data)
    remove('temp_xs.csv')


def test_from_multigroup_flux():
    energies = [0., 6.25e-1, 5.53e3, 8.21e5, 2.e7]
    flux = [1.1e-7, 1.2e-6, 1.3e-5, 1.4e-4]
    chain_file = Path(__file__).parents[1] / 'chain_simple.xml'
    kwargs = {'multigroup_flux': flux, 'chain_file': chain_file}

    # test with energy group structure from string
    microxs = MicroXS.from_multigroup_flux(energies='CASMO-4', **kwargs)
    assert isinstance(microxs, MicroXS)

    # test with energy group structure as floats
    microxs = MicroXS.from_multigroup_flux(energies=energies, **kwargs)
    assert isinstance(microxs, MicroXS)

    # test with nuclides provided
    microxs = MicroXS.from_multigroup_flux(
        energies=energies, nuclides=['Gd157', 'H1'], **kwargs
    )
    assert isinstance(microxs, MicroXS)
    assert microxs.nuclides == ['Gd157', 'H1']

    # test with reactions provided
    microxs = MicroXS.from_multigroup_flux(
        energies=energies, reactions=['fission', '(n,2n)'], **kwargs
    )
    assert isinstance(microxs, MicroXS)
    assert microxs.reactions == ['fission', '(n,2n)']


def test_multigroup_flux_same():
    chain_file = Path(__file__).parents[1] / 'chain_simple.xml'

    # Generate micro XS based on 4-group flux
    energies = [0., 6.25e-1, 5.53e3, 8.21e5, 2.e7]
    flux_per_ev = [0.3, 0.3, 1.0, 1.0]
    flux = flux_per_ev * np.diff(energies)
    flux_sum = flux.sum()
    microxs_4g = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=chain_file)

    # from_multigroup_flux should not modify the flux
    assert flux.sum() == flux_sum

    # Generate micro XS based on 2-group flux, where the boundaries line up with
    # the 4 group flux and have the same flux per eV across the full energy
    # range
    energies = [0., 5.53e3, 2.0e7]
    flux_per_ev = [0.3, 1.0]
    flux = flux_per_ev * np.diff(energies)
    microxs_2g = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=chain_file)

    assert microxs_4g.data == pytest.approx(microxs_2g.data)


def test_microxs_zero_flux():
    chain_file = Path(__file__).parents[1] / 'chain_simple.xml'

    # Generate micro XS based on zero flux
    energies = [0., 6.25e-1, 5.53e3, 8.21e5, 2.e7]
    flux = [0.0, 0.0, 0.0, 0.0]
    microxs = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=chain_file)

    # All microscopic cross sections should be zero
    assert np.all(microxs.data == 0.0)


@pytest.mark.parametrize('include_model_tallies', [None, False, True])
def test_hybrid_tally_setup(include_model_tallies):
    """In hybrid mode a 1-group RR tally is added alongside the flux tally."""
    # Create a simple model with one material and a few nuclides for testing
    model = openmc.Model()
    mat = openmc.Material(components={'U235': 1.0, 'O16': 2.0})
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=mat)
    model.geometry = openmc.Geometry([cell])
    model.settings.batches = 2
    model.settings.particles = 10

    # Define 2-group energy structure for the test
    energies = [0., 0.625, 2.0e7]

    user_tally = openmc.Tally(name='User tally')
    user_tally.scores = ['flux']
    model.tallies = [user_tally]

    # Function to replace Model.run and capture the tallies that were created
    captured = {}
    def capture_run(**kwargs):
        captured['tallies'] = list(model.tallies)
        raise StopIteration

    # Call get_microxs_and_flux but replace Model.run with a function that
    # captures the tallies and raises StopIteration to exit early
    kwargs = {
        'nuclides': ['U235', 'O16'],
        'reactions': ['fission', '(n,gamma)'],
        'energies': energies,
        'reaction_rate_mode': 'flux',
        'reaction_rate_opts': {
            'nuclides': ['U235'], 'reactions': ['fission']
        },
        'chain_file': CHAIN_FILE,
    }
    if include_model_tallies is not None:
        kwargs['include_model_tallies'] = include_model_tallies
    with patch.object(model, 'run', side_effect=capture_run):
        with pytest.raises(StopIteration):
            get_microxs_and_flux(model, [mat], **kwargs)

    # Check that both tallies were created with the expected properties
    tally_names = [t.name for t in captured['tallies']]
    assert ('User tally' in tally_names) is (include_model_tallies is True)
    assert 'MicroXS flux 0' in tally_names
    assert 'MicroXS RR 0' in tally_names

    # Check that the RR tally has the expected nuclides and reactions
    rr = next(t for t in captured['tallies'] if t.name == 'MicroXS RR 0')
    assert rr.nuclides == ['U235']
    assert rr.scores == ['fission']

    # RR tally must use a 1-group energy filter spanning the full energy range
    ef = next(f for f in rr.filters if isinstance(f, openmc.EnergyFilter))
    assert len(ef.values) == 2
    assert ef.values[0] == pytest.approx(energies[0])
    assert ef.values[-1] == pytest.approx(energies[-1])


@pytest.mark.parametrize('split_xml', [False, True])
@pytest.mark.parametrize('mode', [
    'explicit', 'settings-xml-reset', 'model-xml-reset', 'auto', 'shared',
])
def test_microxs_weight_window_mesh(tmp_path, monkeypatch, mode, split_xml):
    monkeypatch.setattr(openmc.MeshBase, 'used_ids', set())
    monkeypatch.setattr(openmc.MeshBase, 'next_id', 1)
    model = openmc.Model(openmc.Geometry([openmc.Cell()]))
    ww_mesh = openmc.RegularMesh()
    ww_mesh.dimension = (1, 1, 1)
    ww_mesh.lower_left = (-1., -1., -1.)
    ww_mesh.upper_right = (1., 1., 1.)
    model.settings.weight_windows = openmc.WeightWindows(
        ww_mesh, [0.5], upper_bound_ratio=2., energy_bounds=[0., 2.e7])

    if mode.endswith('-reset'):
        original = tmp_path / 'original'
        model.export_to_xml(original)
        openmc.MeshBase.reset_ids()
        if mode == 'settings-xml-reset':
            model.settings = openmc.Settings.from_xml(original / 'settings.xml')
        else:
            model = openmc.Model.from_xml(
                geometry=original / 'geometry.xml',
                materials=original / 'materials.xml',
                settings=original / 'settings.xml',
                tallies=original / 'tallies.xml')
        ww_mesh = model.settings.weight_windows[0].mesh
        assert openmc.RegularMesh().id != ww_mesh.id
        openmc.MeshBase.reset_ids()

    assert ww_mesh.id == 1
    if mode == 'shared':
        domain = ww_mesh
    else:
        if mode == 'explicit':
            with pytest.warns(IDWarning, match=r'MeshBase.*id=1'):
                domain = openmc.RegularMesh(mesh_id=1)
        else:
            domain = openmc.RegularMesh()
        assert (domain.id == ww_mesh.id) is (mode != 'auto')
        domain.dimension = (2, 1, 1)
        domain.lower_left = ww_mesh.lower_left
        domain.upper_right = ww_mesh.upper_right

    def export(**kwargs):
        if split_xml:
            model.export_to_xml(tmp_path)
        else:
            model.export_to_model_xml(tmp_path)
        raise StopIteration

    monkeypatch.setattr(model, 'run', export)
    with pytest.raises(StopIteration):
        get_microxs_and_flux(
            model, domain, nuclides=['H1'], reactions=['(n,gamma)'],
            chain_file=CHAIN_FILE)

    if split_xml:
        settings = ET.parse(tmp_path / 'settings.xml').getroot()
        tallies = ET.parse(tmp_path / 'tallies.xml').getroot()
    else:
        root = ET.parse(tmp_path / 'model.xml').getroot()
        settings, tallies = root.find('settings'), root.find('tallies')
    ww_id = settings.findtext('weight_windows/mesh')
    flux_id = tallies.findtext("filter[@type='mesh']/bins")
    meshes = {m.get('id'): m for m in settings.findall('mesh')}
    for mesh in tallies.findall('mesh'):
        meshes.setdefault(mesh.get('id'), mesh)

    assert meshes[ww_id].findtext('dimension') == '1 1 1'
    assert meshes[flux_id].findtext('dimension') == (
        '1 1 1' if mode == 'shared' else '2 1 1')
    assert (flux_id == ww_id) is (mode == 'shared')
    if not split_xml:
        assert len(settings.findall('mesh') + tallies.findall('mesh')) == (
            1 if mode == 'shared' else 2)
    assert model.settings.weight_windows[0].mesh is ww_mesh
    assert ww_mesh.id == 1
    assert domain.id == (2 if mode == 'auto' else 1)


@pytest.mark.parametrize('kind', [
    'rectilinear', 'cylindrical', 'spherical', 'unstructured',
])
@pytest.mark.parametrize('with_weight_windows', [False, True])
def test_microxs_mesh_collision_with_model_tally(
    tmp_path, monkeypatch, kind, with_weight_windows
):
    model = openmc.Model(openmc.Geometry([openmc.Cell()]))
    existing_mesh = openmc.RegularMesh()
    existing_mesh.dimension = (1, 1, 1)
    existing_mesh.lower_left = (-1., -1., -1.)
    existing_mesh.upper_right = (1., 1., 1.)
    existing_filter = openmc.MeshFilter(existing_mesh)
    model.tallies = [openmc.Tally()]
    model.tallies[0].filters = [existing_filter]
    model.tallies[0].scores = ['flux']
    original_tallies = model.tallies
    original_tally = model.tallies[0]

    if with_weight_windows:
        other_mesh = openmc.RegularMesh()
        other_mesh.dimension = (1, 2, 1)
        other_mesh.lower_left = existing_mesh.lower_left
        other_mesh.upper_right = existing_mesh.upper_right
        model.settings.weight_windows = [
            openmc.WeightWindows(
                mesh, [0.5] * mesh.n_elements, upper_bound_ratio=2.,
                energy_bounds=[0., 2.e7])
            for mesh in (other_mesh, existing_mesh)
        ]
    original_windows = model.settings.weight_windows

    mesh_id = existing_mesh.id
    constructors = {
        'rectilinear': (openmc.RectilinearMesh, ()),
        'cylindrical': (openmc.CylindricalMesh, ([0., 1.], [0., 1.5])),
        'spherical': (openmc.SphericalMesh, ([0., 1.],)),
        'unstructured': (openmc.UnstructuredMesh, ('mesh.e',)),
    }
    mesh_type, args = constructors[kind]
    with pytest.warns(IDWarning, match=rf'MeshBase.*id={mesh_id}'):
        domain = mesh_type(*args, mesh_id=mesh_id)
    if kind == 'rectilinear':
        domain.x_grid = domain.y_grid = domain.z_grid = [0., 1.]
    domain_filter = openmc.MeshFilter(domain)
    domain_filter.translation = (1., 2., 3.)
    domain_filter_id = domain_filter.id

    def export(**kwargs):
        model.export_to_model_xml(tmp_path)
        raise StopIteration

    monkeypatch.setattr(model, 'run', export)
    with pytest.raises(StopIteration):
        get_microxs_and_flux(
            model, domain_filter, nuclides=['H1'], reactions=['(n,gamma)'],
            chain_file=CHAIN_FILE, include_model_tallies=True)

    root = ET.parse(tmp_path / 'model.xml').getroot()
    tallies = root.find('tallies')
    flux_tally = tallies.find("tally[@name='MicroXS flux 0']")
    flux_filter = tallies.find(f"filter[@id='{flux_tally.findtext('filters')}']")
    flux_mesh_id = flux_filter.findtext('bins')
    assert flux_mesh_id != str(mesh_id)
    assert tallies.find(f"mesh[@id='{flux_mesh_id}']").get('type') == kind
    assert flux_filter.get('translation') == '1.0 2.0 3.0'
    if with_weight_windows:
        settings = root.find('settings')
        assert [w.findtext('mesh') for w in settings.findall('weight_windows')] == [
            str(other_mesh.id), str(mesh_id)]
        assert settings.find(f"mesh[@id='{mesh_id}']").findtext('dimension') == '1 1 1'
    assert model.settings.weight_windows is original_windows
    assert model.tallies is original_tallies
    assert len(model.tallies) == 1
    assert model.tallies[0] is original_tally
    assert original_tally.filters[0] is existing_filter
    assert domain_filter.id == domain_filter_id
    assert domain_filter.mesh is domain
    assert domain.id == mesh_id


@pytest.mark.parametrize('split_xml', [False, True])
@pytest.mark.parametrize('energies', [None, [0., 1., 2.e7]],
                         ids=['no-energy', 'energy'])
@pytest.mark.parametrize('include_model_tallies', [True, False],
                         ids=['retained-tally', 'domain-sequence'])
def test_microxs_filter_ids_after_reset(
    tmp_path, monkeypatch, include_model_tallies, energies, split_xml
):
    # Isolate the ID pools so exercising the public reset is local to this test.
    for cls in IDManagerMixin.__subclasses__():
        monkeypatch.setattr(cls, 'used_ids', set())
        monkeypatch.setattr(cls, 'next_id', 1)

    cell = openmc.Cell()
    model = openmc.Model(openmc.Geometry([cell]))
    ww_mesh = openmc.RegularMesh(mesh_id=1)
    ww_mesh.dimension = (1, 1, 1)
    ww_mesh.lower_left = (-1., -1., -1.)
    ww_mesh.upper_right = (1., 1., 1.)
    windows = openmc.WeightWindows(
        ww_mesh, [0.5], upper_bound_ratio=2., energy_bounds=[0., 2.e7])
    model.settings.weight_windows = windows

    if include_model_tallies:
        existing_filter = openmc.MeshFilter(ww_mesh, filter_id=2)
        tally = openmc.Tally(tally_id=10)
        tally.filters = [existing_filter]
        tally.scores = ['flux']
        model.tallies = [tally]
    else:
        existing_filter = openmc.CellFilter(cell, filter_id=1)

    original_tallies = model.tallies
    existing_filter_id = existing_filter.id
    openmc.reset_auto_ids()
    assert existing_filter_id not in openmc.Filter.used_ids
    domain = openmc.RegularMesh()
    domain.dimension = (2, 1, 1)
    domain.lower_left = ww_mesh.lower_left
    domain.upper_right = ww_mesh.upper_right
    assert domain.id == ww_mesh.id == 1
    domains = domain
    if not include_model_tallies:
        domain_filter = openmc.MeshFilter(domain, filter_id=50)
        domains = [domain_filter, existing_filter]

    def export(**kwargs):
        if split_xml:
            model.export_to_xml(tmp_path)
        else:
            model.export_to_model_xml(tmp_path)
        raise StopIteration

    monkeypatch.setattr(model, 'run', export)
    with pytest.raises(StopIteration):
        get_microxs_and_flux(
            model, domains, nuclides=['H1'], reactions=['(n,gamma)'],
            energies=energies, chain_file=CHAIN_FILE,
            include_model_tallies=include_model_tallies)

    if split_xml:
        settings = ET.parse(tmp_path / 'settings.xml').getroot()
        tallies = ET.parse(tmp_path / 'tallies.xml').getroot()
    else:
        root = ET.parse(tmp_path / 'model.xml').getroot()
        settings, tallies = root.find('settings'), root.find('tallies')

    filter_elements = tallies.findall('filter')
    filters = {f.get('id'): f for f in filter_elements}
    assert len(filters) == len(filter_elements) == (2 if energies is None else 3)
    for tally_xml in tallies.findall('tally'):
        ids = tally_xml.findtext('filters').split()
        assert len(set(ids)) == len(ids)
        assert all(i in filters for i in ids)
    flux_ids = tallies.findtext("tally[@name='MicroXS flux 0']/filters").split()
    assert filters[flux_ids[0]].get('type') == 'mesh'
    if energies is not None:
        energy_xml = filters[flux_ids[1]]
        assert energy_xml.get('type') == 'energy'
        assert [float(x) for x in energy_xml.findtext('bins').split()] == energies
    assert existing_filter_id != int(flux_ids[0])
    assert int(flux_ids[0]) in openmc.Filter.used_ids
    assert existing_filter_id not in openmc.Filter.used_ids

    ww_id = settings.findtext('weight_windows/mesh')
    assert ww_id == '1'
    assert settings.find(f"mesh[@id='{ww_id}']").findtext('dimension') == '1 1 1'
    mesh_id = filters[flux_ids[0]].findtext('bins')
    assert mesh_id != ww_id
    assert tallies.find(f"mesh[@id='{mesh_id}']").findtext('dimension') == '2 1 1'
    assert model.settings.weight_windows[0] is windows
    assert windows.mesh is ww_mesh
    assert ww_mesh.id == domain.id == 1
    assert model.tallies is original_tallies
    assert existing_filter.id == existing_filter_id
    if include_model_tallies:
        assert model.tallies[0].filters[0] is existing_filter
        assert existing_filter.mesh is ww_mesh
        assert tallies.findtext("tally[@id='10']/filters") == str(existing_filter_id)
        assert filters[str(existing_filter_id)].findtext('bins') == ww_id
    else:
        assert not model.tallies
        assert domain_filter.id == 50
        assert domain_filter.id in openmc.Filter.used_ids
        assert domain_filter.mesh is domain
        cell_flux_ids = tallies.findtext("tally[@name='MicroXS flux 1']/filters")
        assert cell_flux_ids.split()[0] == str(existing_filter_id)
        assert filters[str(existing_filter_id)].findtext('bins') == str(cell.id)


@pytest.mark.parametrize('split_xml', [False, True], ids=['combined', 'split'])
@pytest.mark.parametrize('energies', [None, [0., 1.e5, 2.e7]],
                         ids=['no-energy', 'energy'])
@pytest.mark.parametrize('retained_id', [1, 2], ids=['flux-id', 'rr-id'])
def test_microxs_tally_ids_after_reset(
    tmp_path, monkeypatch, retained_id, energies, split_xml
):
    for cls in IDManagerMixin.__subclasses__():
        monkeypatch.setattr(cls, 'used_ids', set())
        monkeypatch.setattr(cls, 'next_id', 1)

    cell = openmc.Cell()
    model = openmc.Model(openmc.Geometry([cell]))
    if retained_id == 2:
        openmc.Tally()  # Consume the ID normally assigned to the flux tally.
    retained_tally = openmc.Tally(name='retained')
    retained_filter = openmc.CellFilter(cell)
    retained_tally.filters = [retained_filter]
    retained_tally.scores = ['flux']
    model.tallies = [retained_tally]
    original_tallies = model.tallies
    assert retained_tally.id == retained_id
    openmc.reset_auto_ids()
    assert retained_id not in openmc.Tally.used_ids

    def export(**kwargs):
        if split_xml:
            model.export_to_xml(tmp_path)
        else:
            model.export_to_model_xml(tmp_path)
        raise StopIteration

    monkeypatch.setattr(model, 'run', export)
    with pytest.raises(StopIteration):
        get_microxs_and_flux(
            model, [cell], nuclides=['H1'], reactions=['(n,gamma)'],
            energies=energies, chain_file=CHAIN_FILE, include_model_tallies=True)

    assert model.tallies is original_tallies
    assert len(model.tallies) == 1
    assert model.tallies[0] is retained_tally
    assert retained_tally.id == retained_id
    assert retained_tally.filters[0] is retained_filter
    assert retained_tally.scores == ['flux']
    if split_xml:
        tallies = ET.parse(tmp_path / 'tallies.xml').getroot()
    else:
        tallies = ET.parse(tmp_path / 'model.xml').getroot().find('tallies')

    entries = tallies.findall('tally')
    ids = [int(entry.get('id')) for entry in entries]
    assert len(set(ids)) == len(entries) == 3, f'Duplicate tally IDs: {ids}'
    retained_xml, flux_xml, rr_xml = entries
    assert int(retained_xml.get('id')) == retained_id
    assert retained_xml.findtext('filters') == str(retained_filter.id)
    assert flux_xml.findtext('scores') == 'flux'
    assert rr_xml.findtext('scores') == '(n,gamma)'
    filter_ids = flux_xml.findtext('filters').split()
    assert filter_ids == rr_xml.findtext('filters').split()
    assert len(filter_ids) == (1 if energies is None else 2)
    assert set(filter_ids) <= {f.get('id') for f in tallies.findall('filter')}
    assert {int(flux_xml.get('id')), int(rr_xml.get('id'))} <= openmc.Tally.used_ids
    assert retained_id not in openmc.Tally.used_ids


def _simple_model():
    model = openmc.Model()
    mat = openmc.Material(components={'H1': 1.0, 'H2': 1.0},
                          density=5.0, density_units='g/cm3')
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=mat)
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 100
    model.settings.batches = 5
    model.settings.run_mode = 'fixed source'
    return model, mat


def test_hybrid_tally_defaults_to_all_nuclides(run_in_tmpdir):
    energies = [0., 0.625, 2.0e7]
    kwargs = {
        'nuclides': ['H1', 'H2'],
        'reactions': ['(n,2n)', '(n,gamma)'],
        'energies': energies,
        'reaction_rate_mode': 'flux',
        'chain_file': CHAIN_FILE,
    }

    model, mat = _simple_model()
    default_fluxes, default_micros = get_microxs_and_flux(
        model, [mat], reaction_rate_opts={'reactions': ['(n,2n)']}, **kwargs
    )

    model, mat = _simple_model()
    explicit_fluxes, explicit_micros = get_microxs_and_flux(
        model, [mat],
        reaction_rate_opts={
            'nuclides': ['H1', 'H2'],
            'reactions': ['(n,2n)']
        },
        **kwargs
    )

    np.testing.assert_allclose(default_fluxes[0], explicit_fluxes[0])
    np.testing.assert_allclose(default_micros[0].data, explicit_micros[0].data)
    assert default_micros[0].nuclides == explicit_micros[0].nuclides
    assert default_micros[0].reactions == explicit_micros[0].reactions


def test_flux_mode_returns_one_group_flux(run_in_tmpdir):
    model, mat = _simple_model()
    fluxes, micros = get_microxs_and_flux(
        model, [mat],
        nuclides=['H1'],
        reactions=['(n,2n)'],
        energies=[0., 0.625, 2.0e7],
        reaction_rate_mode='flux',
        chain_file=CHAIN_FILE,
    )

    assert fluxes[0].shape == (1,)
    assert micros[0].data.shape == (1, 1, 1)
    assert fluxes[0][0] > 0.0

# ---------------------------------------------------------------------------
# Tests for MicroXS.merge()
# ---------------------------------------------------------------------------

def _make_microxs(nuclides, reactions, values, groups=1):
    """Helper: build a MicroXS from a flat list of values (nuclide-major order)."""
    data = np.array(values, dtype=float).reshape(
        len(nuclides), len(reactions), groups)
    return MicroXS(data, nuclides, reactions)


def test_merge_disjoint():
    """Merging two MicroXS with no overlapping nuclides or reactions."""
    m1 = _make_microxs(['U235', 'U238'], ['fission', '(n,gamma)'],
                        [1., 2., 3., 4.])
    m2 = _make_microxs(['Pu239'], ['(n,2n)'], [5.])

    merged = m1.merge(m2)

    assert merged.nuclides == ['U235', 'U238', 'Pu239']
    assert merged.reactions == ['fission', '(n,gamma)', '(n,2n)']
    assert merged.data.shape == (3, 3, 1)

    # Self data preserved
    assert merged['U235', 'fission'] == pytest.approx([1.])
    assert merged['U238', '(n,gamma)'] == pytest.approx([4.])
    # New data from other
    assert merged['Pu239', '(n,2n)'] == pytest.approx([5.])
    # Cross-terms that had no data should be zero
    assert merged['U235', '(n,2n)'] == pytest.approx([0.])
    assert merged['Pu239', 'fission'] == pytest.approx([0.])


def test_merge_prefer_other():
    """prefer='other': other overwrites shared entries, adds new ones."""
    # m1: U235 and U238, reactions fission and (n,gamma)
    m1 = _make_microxs(['U235', 'U238'], ['fission', '(n,gamma)'],
                        [1., 2., 3., 4.])
    # m2: only U235, reactions fission (overlap) and (n,2n) (new)
    m2 = _make_microxs(['U235'], ['fission', '(n,2n)'], [9., 5.])

    merged = m1.merge(m2)

    # Nuclide/reaction sets
    assert set(merged.nuclides) == {'U235', 'U238'}
    assert set(merged.reactions) == {'fission', '(n,gamma)', '(n,2n)'}

    # U235/fission overwritten by m2
    assert merged['U235', 'fission'] == pytest.approx([9.])
    # U235/(n,gamma) untouched
    assert merged['U235', '(n,gamma)'] == pytest.approx([2.])
    # New reaction added from m2
    assert merged['U235', '(n,2n)'] == pytest.approx([5.])
    # U238 data preserved
    assert merged['U238', 'fission'] == pytest.approx([3.])
    assert merged['U238', '(n,gamma)'] == pytest.approx([4.])
    # U238/(n,2n) not in either → zero
    assert merged['U238', '(n,2n)'] == pytest.approx([0.])


def test_merge_prefer_self():
    """prefer='self': shared pairs keep self's value; new entries use other's value."""
    m1 = _make_microxs(['U235', 'U238'], ['fission', '(n,gamma)'],
                        [1., 2., 3., 4.])
    m2 = _make_microxs(['U235'], ['fission', '(n,2n)'], [9., 5.])

    merged = m1.merge(m2, prefer='self')

    # U235/fission: self wins
    assert merged['U235', 'fission'] == pytest.approx([1.])
    # U235/(n,2n): other used
    assert merged['U235', '(n,2n)'] == pytest.approx([5.])
