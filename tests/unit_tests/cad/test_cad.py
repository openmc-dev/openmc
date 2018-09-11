
from subprocess import call

import os
import sys


import pytest
pytestmark = pytest.mark.skipif(
    os.environ.get('CAD') != 'y',
    reason="CAD build is not enabled.")

def test_cad():
    # make sure cwd is in the python system path
    d = os.path.abspath(os.path.dirname(__file__))
    # move here
    os.chdir(d)
    # run test
    assert not call(['openmc','-s','1'])
    
