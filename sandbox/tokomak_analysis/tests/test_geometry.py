import numpy as np
import pytest
from dataclasses import replace
from sandbox.tokomak_analysis import Config,build_model
from sandbox.tokomak_analysis.validation import walk
from sandbox.tokomak_analysis.geometry.primitives import rz_halfspace,polygon_volume


def test_signed_rz_halfplanes():
    rng=np.random.default_rng(9)
    for angle in np.linspace(0,2*np.pi,17):
        normal=np.array([np.cos(angle),np.sin(angle)]);constant=normal[0]*910
        region=rz_halfspace(normal,constant)
        for r,z in rng.uniform([1,-1000],[1800,1000],(20,2)):
            assert ((r,0,z) in region)==(normal@[r,z]>=constant)


def test_random_hierarchy_coverage(assembly):
    rng=np.random.default_rng(87)
    for p in rng.uniform([-1800,-1800,-1150],[1800,1800,1150],(2000,3)):
        if np.linalg.norm(p)<assembly.config.world_radius:assert walk(assembly.model.geometry.root_universe,p)


def test_exact_plasma_volume(assembly):
    volume=polygon_volume(assembly.contours['plasma'].points)
    assert 2000e6<volume<3000e6


def test_native_csg_only_and_xml(assembly,tmp_path):
    assembly.model.export_to_xml(directory=tmp_path)
    assert (tmp_path/'geometry.xml').exists()
    assert all(s.type in ('plane','x-plane','y-plane','z-plane','sphere','z-cylinder','z-cone') for s in assembly.model.geometry.get_all_surfaces().values())
    assert not any(cell.fill_type=='lattice' for cell in assembly.model.geometry.get_all_cells().values())


def test_sector_periodicity(assembly):
    angle=2*np.pi/assembly.config.tf_coils;p=np.array([1300.,20.,100.])
    roles=[]
    for j in range(assembly.config.tf_coils):
        t=j*angle;q=np.array([p[0]*np.cos(t)-p[1]*np.sin(t),p[0]*np.sin(t)+p[1]*np.cos(t),p[2]])
        leaf=walk(assembly.model.geometry.root_universe,q)[-1]
        roles.append(assembly.registry.records[leaf.id]['component'])
    assert len(set(roles))==1


@pytest.mark.parametrize('settings',[dict(fidelity=4),dict(breeder_packing=1.2),dict(enrichment=(.6,)),dict(poloidal_modules=3)])
def test_invalid_configuration(settings):
    with pytest.raises(ValueError):replace(Config(),**settings).validate()


def test_validate_cli_dispatch(monkeypatch, tmp_path, assembly):
    """The standalone validation command must not require transport locals."""
    import sys
    import sandbox.tokomak_analysis.main as cli
    import sandbox.tokomak_analysis.validation as validation
    called = []
    monkeypatch.setattr(cli, 'build_model', lambda *args, **kwargs: assembly)
    monkeypatch.setattr(cli, 'export_inputs', lambda *args: None)
    monkeypatch.setattr(validation, 'sample_geometry',
                        lambda model, path, count: called.append((path, count)))
    monkeypatch.setattr(sys, 'argv', ['main', '--validate', '--output', str(tmp_path)])
    cli.main()
    assert called == [(tmp_path, 20000)]
