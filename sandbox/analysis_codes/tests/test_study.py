"""Study separation, scoring, regression and provenance checks."""
from pathlib import Path
import ast
import json
import shutil
import numpy as np
import pytest
import openmc
from sandbox.analysis_codes import context,pce,sobol
from sandbox.analysis_codes.analysis import paired_morris_errors
from sandbox.analysis_codes.workflow import guard_request

ROOT=Path(__file__).resolve().parents[2]


def test_no_external_analysis_imports():
    forbidden=('Li_uq','FLiBe_uq','FLiBe_analysis','nuclear_data_sensitivity','mt205_investigation')
    for file in (ROOT/'analysis_codes').glob('*.py'):
        for node in ast.walk(ast.parse(file.read_text())):
            if isinstance(node,ast.ImportFrom):
                assert not any((node.module or '').startswith(name) for name in forbidden),file
            if isinstance(node,ast.Import):
                assert not any(n.name.startswith(forbidden) for n in node.names),file


def test_response_weights():
    context.load_study(ROOT/'studies'/'FLiBe_analysis')
    openmc.reset_auto_ids();model=context.build_model();rsp=context.responses(model)
    assert set(rsp)=={'total','Li6','Li7','Be9','F19'}
    tally=model.tallies[-1]
    assert '205' not in tally.scores
    weights=np.asarray(rsp['Be9'].score_weights)
    assert weights[:,tally.scores.index('105')].sum()==0
    for mt in (700,701):assert weights[:,tally.scores.index(str(mt))].sum()==1
    np.testing.assert_array_equal(rsp['total'].score_weights,
        sum(np.asarray(rsp[n].score_weights) for n in rsp if n!='total'))


def test_study_load_and_guard(tmp_path):
    study=tmp_path/'A study';study.mkdir()
    for file in ('model.py','config.py'):shutil.copy2(ROOT/'studies'/'FLiBe_analysis'/file,study/file)
    loaded=context.load_study(study)
    guard_request(loaded)
    previous=(loaded.output/'analysis_request.json').read_bytes()
    loaded.config.PCE_BACKEND='openturns'
    guard_request(loaded)  # Statistical backend does not change transport inputs.
    loaded.config.SCREEN_PARTICLES+=1
    with pytest.raises(ValueError,match='new --output'):guard_request(loaded)
    assert (loaded.output/'analysis_request.json').read_bytes()==previous
    other=tmp_path/'another';other.mkdir()
    (other/'config.py').write_text('VALUE=123\n')
    (other/'model.py').write_text('from .config import VALUE\ndef build_model(particles,batches): return VALUE\n')
    context.load_study(other)
    assert context.build_model()==123
    assert loaded.config.SCREEN_PARTICLES==1001


def test_paired_morris():
    from types import SimpleNamespace
    b=np.zeros((1,2,10,1));b[0,0,:,0]=np.arange(10);b[0,1,:,0]=np.arange(10)+.3
    design=SimpleNamespace(changed=np.array([[0]]),steps=np.array([[.6]]))
    assert paired_morris_errors(design,b)[0,0,0]<1e-14


def test_orthogonality_and_sobol():
    x,w=np.polynomial.legendre.leggauss(8)
    u=np.array(np.meshgrid((x+1)/2,(x+1)/2)).reshape(2,-1).T
    powers=pce.indices(2,2);a=pce.matrix(u,powers)
    np.testing.assert_allclose(a.T@(np.outer(w,w).ravel()[:,None]/4*a),np.eye(a.shape[1]),atol=1e-13)
    x=np.random.default_rng(7).random((120,2));powers=pce.indices(2,2)
    c=np.zeros(len(powers))
    for mode,value in [((0,0),3),((1,0),2),((0,1),1),((1,1),3)]:
        c[np.where((powers==mode).all(axis=1))[0][0]]=value
    model=pce.fit(x,pce.matrix(x,powers)@c,np.full(120,1e-6),order=2)
    result=sobol.analyse(model,groups={'both':[0,1]})
    np.testing.assert_allclose(result['first'][:,0],[4/14,1/14],atol=1e-12)
    np.testing.assert_allclose(result['total'][:,0],[13/14,10/14],atol=1e-12)
    assert result['groups']['both']['total'][0]==pytest.approx(1)
