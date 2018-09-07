
import pytest

from subprocess import call

from tests.unit_tests import assert_false

def test_cad():
    assert not call(["openmc"])
    
