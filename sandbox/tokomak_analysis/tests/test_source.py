from dataclasses import replace
import numpy as np
from sandbox.tokomak_analysis.config import Config
from sandbox.tokomak_analysis.neutronics.source import sample_profile,radial_coordinate,make_source


def test_source_power_and_profile(assembly):
    c=Config();assert abs(c.neutron_rate-7.092623948e20)/c.neutron_rate<1e-9
    boundary=assembly.contours['plasma'].points
    uniform=sample_profile(c,boundary,3000,123)
    peaked=sample_profile(replace(c,radial_profile_power=2),boundary,3000,123)
    radius=lambda p:radial_coordinate(np.column_stack([np.hypot(p[:,0],p[:,1]),p[:,2]]),boundary,c.major_radius)
    assert (radius(uniform)<=1+1e-12).all()
    assert radius(peaked).mean()<radius(uniform).mean()


def test_native_profile_source_exists(assembly):
    c=replace(Config(),radial_profile_power=1)
    src=make_source(c,next(cell for cell in assembly.registry.cells.values() if cell.name=='Plasma source domain'),assembly.contours['plasma'].points)
    assert src.to_xml_element().get('type')=='tokamak'


def test_native_source_samples_stay_in_plasma(tmp_path,monkeypatch):
    from sandbox.tokomak_analysis import build_model
    from sandbox.tokomak_analysis.validation import walk
    monkeypatch.chdir(tmp_path)
    for power in (0.,1.):
        c=replace(Config(),radial_profile_power=power,enable_mesh=False)
        a=build_model(c,tallies=False)
        sites=a.model.sample_external_source(500,prn_seed=765)
        assert all(a.registry.records[walk(a.model.geometry.root_universe,site.r)[-1].id]['component']=='plasma' for site in sites)


def test_optional_poloidal_source_bank(tmp_path,monkeypatch):
    from sandbox.tokomak_analysis import build_model
    from sandbox.tokomak_analysis.validation import walk
    monkeypatch.chdir(tmp_path)
    c=replace(Config(),poloidal_asymmetry=.4,radial_profile_power=1.,source_bank_size=1000,enable_mesh=False)
    a=build_model(c,tallies=False,output=tmp_path)
    assert (tmp_path/'plasma_source.h5').exists()
    sites=a.model.sample_external_source(200,prn_seed=765)
    assert all(a.registry.records[walk(a.model.geometry.root_universe,site.r)[-1].id]['component']=='plasma' for site in sites)
    np.testing.assert_allclose([np.linalg.norm(s.u) for s in sites],1,atol=1e-12)
