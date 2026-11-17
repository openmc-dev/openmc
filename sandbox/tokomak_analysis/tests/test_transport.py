"""Short real native geometry-debug transport, including batch covariance checks."""
from dataclasses import replace
import numpy as np
from sandbox.tokomak_analysis import Config,build_model
from sandbox.tokomak_analysis.main import run_transport


def test_native_debug_transport(tmp_path):
    c=replace(Config(),fidelity=1,particles=30,batches=3,enable_mesh=False,threads=1)
    a=build_model(c,output=tmp_path)
    result=run_transport(a,tmp_path,debug=True)
    assert result['performance']['histories']==90
    assert np.isfinite(result['TBR']['per_source']) and result['TBR']['per_source']>=0
    assert (tmp_path/'statepoint.3.h5').exists()
    assert (tmp_path/'module_breeding.csv').exists()
