
from subprocess import call

def test_cad():
    assert not call(['openmc'])
    
