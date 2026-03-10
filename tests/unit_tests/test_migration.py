import openmc
import pytest
import numpy as np

from tests import cdtemp

def sphere_model(radius, boundary_type):
    
    model = openmc.Model()

    # Material
    material = openmc.Material(name="Hydrogen")
    material.add_nuclide("H1", 1.0)
    material.set_density('g/cm3', 1.0)
    material.add_s_alpha_beta('c_H_in_H2O')

    # Geometry
    sphere = openmc.Sphere(r=radius, boundary_type=boundary_type)
    cell = openmc.Cell(region=-sphere, fill=material)
    model.geometry = openmc.Geometry([cell])

    # Settings
    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource()


    # Tally
    tally = openmc.Tally()
    tally.scores = ["migration-area"]
    model.tallies = [tally]

    return model, tally

def test_reflection_is_equivalent_to_large_model(run_in_tmpdir):
    openmc.reset_auto_ids()
    refl_model, tally1 = sphere_model(2.5, "reflective")
    large_model, tally2 = sphere_model(100, "vacuum")

    with cdtemp():
        refl_model.run(apply_tally_results=True)
        mean1, std1 = (tally1.mean.squeeze(),
                       tally1.std_dev.squeeze())
    
    with cdtemp():
        large_model.run(apply_tally_results=True)
        mean2, std2 = (tally2.mean.squeeze(),
                       tally2.std_dev.squeeze())
    std = (std1**2+std2**2)**0.5
    diff = np.abs(mean2-mean1)
        
    assert diff <= 3*std
