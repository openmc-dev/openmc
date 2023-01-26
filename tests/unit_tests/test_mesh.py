import openmc
import pytest

def test_raises_error_when_flat():
    with pytest.raises(ValueError):
        mesh = openmc.RegularMesh()
        mesh.dimension = [50, 50, 1]
        mesh.lower_left = [-25, -25, 0]
        mesh.upper_right = [25, 25, 0]