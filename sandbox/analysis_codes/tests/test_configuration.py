"""Configured bounds must reach native candidates and exported transport data."""
from dataclasses import asdict
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
from sandbox.analysis_codes.registry import Parameter
from sandbox.analysis_codes.perturbations import shape
from sandbox.analysis_codes.family_workflow import high_transport
from openmc.lib.perturbation import _areas


@pytest.mark.parametrize('limit',[.05,.10,.35,.8])
@pytest.mark.parametrize('interpolation',[1,2])
def test_pdf_bounds_and_normalization(limit,interpolation):
    table=dict(x=np.array([-1.,-.5,0.,.5,1.]),pdf=np.full(5,.5),
               cdf=np.linspace(0,1,5),interpolation=interpolation)
    law=dict(type='tabular-angle',incident_energy=np.array([2e6]),tables=[table])
    params=[Parameter('Li7',2,'angle',mode=i,limit=limit) for i in (1,2)]
    for values in ([limit,limit],[-limit,-limit],[limit,0]):
        changed=shape(law,params,np.asarray(values))['tables'][0]
        assert _areas(changed['x'],changed['pdf'],interpolation).sum()==pytest.approx(1)
        assert np.max(abs(changed['pdf']/.5-1))<=limit+1e-12
        assert changed['cdf'][-1]==pytest.approx(1)


def test_high_transport_uses_registry_limits_and_run_settings(tmp_path,monkeypatch):
    from sandbox.analysis_codes import family_workflow as module
    params=[Parameter('Li6',105,'xs',limit=.10),Parameter('Li7',2,'angle',mode=1,limit=.35)]
    registry=[dict(**asdict(p),parameter_name=p.identity) for p in params]
    selected=pd.DataFrame({'parameter':[p.identity for p in params]})
    seen=[];built=[]
    class FakeRunner:
        def __init__(self,spec,output,cap):
            assert cap==12345
            self.provenance={};self.history_budget=SimpleNamespace(used=lambda:0)
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def bind(self,parameters):assert parameters==params
        def evaluate(self,values,*,stage,seed):
            seen.append((np.asarray(values),seed))
            return {'responses':{'total':([1.],[.01])}}
    def factory(particles,batches):
        built.append((particles,batches));return object()
    monkeypatch.setattr(module,'BudgetedRunner',FakeRunner)
    high_transport(tmp_path,selected,registry,training=7,validation=3,particles=41,batches=6,
        seed=73,history_cap=12345,model_factory=factory,response_factory=lambda model:{'total':object()},
        spec_factory=lambda model,r:SimpleNamespace(responses=r))
    assert built==[(41,6)]
    bundle=np.load(tmp_path/'high/analysis_arrays.npz')
    expected=np.vstack([np.zeros(2),(2*bundle['training']-1)*[.1,.35],(2*bundle['validation']-1)*[.1,.35]])
    np.testing.assert_allclose([v for v,seed in seen],expected)
    data=pd.read_csv(tmp_path/'high/pce_training.csv')
    np.testing.assert_allclose(data[[p.identity for p in params]],expected[1:8])
    assert (data.histories==246).all()
    assert data.seed.tolist()==list(range(300073,300080))


def test_configured_output_directory(tmp_path):
    from sandbox.analysis_codes.context import load_study
    (tmp_path/'config.py').write_text("RESULTS_DIRECTORY='my_results'\n")
    (tmp_path/'model.py').write_text('def build_model(particles,batches): return None\n')
    study=load_study(tmp_path)
    assert study.output==tmp_path/'my_results'
    study=load_study(tmp_path,tmp_path/'override')
    assert study.output==tmp_path/'override'


@pytest.mark.parametrize('study_name',['Li_analysis','FLiBe_analysis'])
def test_editable_model_settings_reach_openmc(study_name):
    from pathlib import Path
    import openmc
    from sandbox.analysis_codes import context
    study=context.load_study(Path(__file__).resolve().parents[2]/'studies'/study_name)
    c=study.config
    c.DENSITY_G_CM3=.71;c.TEMPERATURE_K=600.;c.LI6_ATOM_FRACTION=.30
    c.RADIUS_CM=15.;c.WORLD_MARGIN_CM=2.
    c.SOURCE_POSITION_CM=(1.,0.,-15.2);c.SOURCE_DIRECTION=(0.,1.,0.)
    c.SOURCE_ENERGY_EV=5e6
    openmc.reset_auto_ids();model=context.build_model(91,7)
    material,=model.geometry.get_all_materials().values()
    assert material.density==.71 and material.temperature==600
    atoms=material.get_nuclide_atom_densities()
    assert atoms['Li6']/(atoms['Li6']+atoms['Li7'])==pytest.approx(.30)
    assert model.settings.particles==91 and model.settings.batches==7
    source=model.settings.source[0]
    np.testing.assert_array_equal(source.space.xyz,[1.,0.,-15.2])
    np.testing.assert_array_equal(source.angle.reference_uvw,[0.,1.,0.])
    np.testing.assert_array_equal(source.energy.x,[5e6])
    assert max(s.r for s in model.geometry.get_all_surfaces().values())==17


def test_morris_confidence_and_resolution_controls():
    from sandbox.analysis_codes.morris import statistics,convergence
    rng=np.random.default_rng(21)
    ee=1+rng.normal(0,.2,(20,2,1));se=np.full_like(ee,.01)
    wide=statistics(ee,se,bootstraps=100,seed=3,confidence_level=.99)
    narrow=statistics(ee,se,bootstraps=100,seed=3,confidence_level=.50)
    assert np.all(wide['mu_star_ci_high']-wide['mu_star_ci_low'] >= narrow['mu_star_ci_high']-narrow['mu_star_ci_low'])
    assert not statistics(ee,se,bootstraps=40,noise_multiplier=1e6)['resolved'].any()
    assert not convergence(ee,se,minimum=21,bootstraps=40)[-1]['stopping_criterion_met']


def test_validation_acceptance_controls():
    from sandbox.analysis_codes.pce import PCE
    from sandbox.analysis_codes.validation import validate
    x=np.linspace(.1,.9,10)[:,None]
    model=PCE(np.array([[0],[1]]),np.array([[2.],[.3]]),np.zeros((1,2,2)),{})
    observed=model.predict(x)+.02
    default=validate(model,x,observed,np.full_like(observed,.01))
    strict=validate(model,x,observed,np.full_like(observed,.01),min_r2=1.,max_noise_mse=0.)
    assert default['accepted'][0] and not strict['accepted'][0]
