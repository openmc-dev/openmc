from pathlib import Path

import pytest
import openmc
from openmc.deplete import Chain, R2SManager


@pytest.fixture
def simple_model_and_mesh(tmp_path):
    # Define two materials: water and Ni
    h2o = openmc.Material()
    h2o.add_nuclide("H1", 2.0)
    h2o.add_nuclide("O16", 1.0)
    h2o.set_density("g/cm3", 1.0)
    nickel = openmc.Material()
    nickel.add_element("Ni", 1.0)
    nickel.set_density("g/cm3", 4.0)

    # Geometry: two half-spaces split by x=0 plane
    left = openmc.XPlane(0.0)
    x_min = openmc.XPlane(-10.0, boundary_type='vacuum')
    x_max = openmc.XPlane(10.0, boundary_type='vacuum')
    y_min = openmc.YPlane(-10.0, boundary_type='vacuum')
    y_max = openmc.YPlane(10.0, boundary_type='vacuum')
    z_min = openmc.ZPlane(-10.0, boundary_type='vacuum')
    z_max = openmc.ZPlane(10.0, boundary_type='vacuum')

    c1 = openmc.Cell(fill=h2o, region=+x_min & -left & +y_min & -y_max & +z_min & -z_max)
    c2 = openmc.Cell(fill=nickel, region=+left & -x_max & +y_min & -y_max & +z_min & -z_max)
    c1.volume = 4000.0
    c2.volume = 4000.0
    geometry = openmc.Geometry([c1, c2])

    # Simple settings with a point source
    settings = openmc.Settings()
    settings.batches = 10
    settings.particles = 1000
    settings.run_mode = 'fixed source'
    settings.source = openmc.IndependentSource()
    model = openmc.Model(geometry, settings=settings)

    mesh = openmc.RegularMesh()
    mesh.lower_left = (-10.0, -10.0, -10.0)
    mesh.upper_right = (10.0, 10.0, 10.0)
    mesh.dimension = (1, 1, 1)
    return model, (c1, c2), mesh


def test_r2s_mesh_expected_output(simple_model_and_mesh, tmp_path):
    model, (c1, c2), mesh = simple_model_and_mesh

    # Use mesh-based domains
    r2s = R2SManager(model, mesh)

    # Use custom reduced chain file for Ni
    chain = Chain.from_xml(Path(__file__).parents[1] / "chain_ni.xml")

    # Run R2S calculation
    outdir = r2s.run(
        timesteps=[(1.0, 'd')],
        source_rates=[1.0],
        photon_time_indices=[1],
        output_dir=tmp_path,
        chain_file=chain,
    )

    # Check directories and files exist
    nt = Path(outdir) / 'neutron_transport'
    assert (nt / 'fluxes.npy').exists()
    assert (nt / 'micros.h5').exists()
    assert (nt / 'mesh_material_volumes.npz').exists()
    act = Path(outdir) / 'activation'
    assert (act / 'depletion_results.h5').exists()
    pt = Path(outdir) / 'photon_transport'
    assert (pt / 'tally_ids.json').exists()
    assert (pt / 'time_1' / 'statepoint.10.h5').exists()

    # Basic results structure checks
    assert len(r2s.results['fluxes']) == 2
    assert len(r2s.results['micros']) == 2
    assert len(r2s.results['mesh_material_volumes']) == 2
    assert len(r2s.results['activation_materials']) == 2
    assert len(r2s.results['depletion_results']) == 2

    # Check activation materials
    amats = r2s.results['activation_materials']
    assert all(m.depletable for m in amats)
    # Volumes preserved
    assert {m.volume for m in amats} == {c1.volume, c2.volume}

    # Check loading results
    r2s_loaded = R2SManager(model, mesh)
    r2s_loaded.load_results(outdir)
    assert len(r2s_loaded.results['fluxes']) == 2
    assert len(r2s_loaded.results['micros']) == 2
    assert len(r2s_loaded.results['mesh_material_volumes']) == 2
    assert len(r2s_loaded.results['activation_materials']) == 2
    assert len(r2s_loaded.results['depletion_results']) == 2


def test_r2s_cell_expected_output(simple_model_and_mesh, tmp_path):
    model, (c1, c2), _ = simple_model_and_mesh

    # Use cell-based domains
    r2s = R2SManager(model, [c1, c2])

    # Use custom reduced chain file for Ni
    chain = Chain.from_xml(Path(__file__).parents[1] / "chain_ni.xml")

    # Run R2S calculation
    bounding_boxes = {c1.id: c1.bounding_box, c2.id: c2.bounding_box}
    outdir = r2s.run(
        timesteps=[(1.0, 'd')],
        source_rates=[1.0],
        photon_time_indices=[1],
        output_dir=tmp_path,
        bounding_boxes=bounding_boxes,
        chain_file=chain
    )

    # Check directories and files exist
    nt = Path(outdir) / 'neutron_transport'
    assert (nt / 'fluxes.npy').exists()
    assert (nt / 'micros.h5').exists()
    act = Path(outdir) / 'activation'
    assert (act / 'depletion_results.h5').exists()
    pt = Path(outdir) / 'photon_transport'
    assert (pt / 'tally_ids.json').exists()
    assert (pt / 'time_1' / 'statepoint.10.h5').exists()

    # Basic results structure checks
    assert len(r2s.results['fluxes']) == 2
    assert len(r2s.results['micros']) == 2
    assert len(r2s.results['activation_materials']) == 2
    assert len(r2s.results['depletion_results']) == 2

    # Check activation materials
    amats = r2s.results['activation_materials']
    assert all(m.depletable for m in amats)
    # Names include cell IDs
    assert any(f"Cell {c1.id}" in m.name for m in amats)
    assert any(f"Cell {c2.id}" in m.name for m in amats)
    # Volumes preserved
    assert {m.volume for m in amats} == {c1.volume, c2.volume}

    # Check loading results
    r2s_loaded = R2SManager(model, [c1, c2])
    r2s_loaded.load_results(outdir)
    assert len(r2s_loaded.results['fluxes']) == 2
    assert len(r2s_loaded.results['micros']) == 2
    assert len(r2s_loaded.results['activation_materials']) == 2
    assert len(r2s_loaded.results['depletion_results']) == 2
