"""Undefined noise diagnostics must not warn or validate constant responses."""
import warnings
import numpy as np
from sandbox.analysis_codes.pce import PCE
from sandbox.analysis_codes.validation import validate


def test_zero_channel_and_resolved_channel():
    x=np.linspace(.1,.9,10)[:,None]
    model=PCE(np.array([[0],[1]]),np.array([[2.,0.],[.3,0.]]),
              np.zeros((2,2,2)),{})
    y=model.predict(x)
    errors=np.column_stack([np.full(10,.01),np.zeros(10)])
    with warnings.catch_warnings():
        warnings.simplefilter('error',RuntimeWarning)
        result=validate(model,x,y,errors)
    np.testing.assert_array_equal(result['noise_diagnostic_samples'],[10,0])
    np.testing.assert_array_equal(result['accepted'],[True,False])
    assert result['noise_mse'][0]==0
    assert np.isnan(result['noise_mse'][1])
    np.testing.assert_array_equal(result['constant_observed_response'],[False,True])
