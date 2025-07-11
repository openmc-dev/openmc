"""This model is based on a numerical example at doi:10.1016/s0306-4549(98)00048-6
"""
import openmc
from openmc.stats import spherical_uniform, Point, Mixture, Watt
import numpy as np
import pytest
import os
import glob
    
from tests.testing_harness import PyAPITestHarness


class SourceTestHarness(PyAPITestHarness):
    def _get_results(self, hash_output=False):
        outstr = super()._get_results(hash_output=hash_output)
        # Read the statepoint file.
        statepoint = glob.glob(self._sp_name)[0]
        with openmc.StatePoint(statepoint) as sp:
            # Write out multiplication.
            outstr += 'multiplication:\n'
            form = '{0:12.6E} {1:12.6E}\n'
            M = sp.multiplication
            outstr += form.format(M.n, M.s)              
            
            # Write out k1.
            outstr += 'k1:\n'
            form = '{0:12.6E} {1:12.6E}\n'
            k1 = sp.k_generation[0]
            outstr += form.format(k1.n, k1.s)   
            
            # Write out kf.
            outstr += 'kf:\n'
            form = '{0:12.6E} {1:12.6E}\n'
            kf = 1-k1/(M-1) 
            outstr += form.format(kf.n, kf.s)
            
            # Write out g*.
            outstr += 'g*:\n'
            form = '{0:12.6E} {1:12.6E}\n'
            g_star = sp.source_efficiency 
            outstr += form.format(g_star.n, g_star.s)
        return outstr   

@pytest.fixture()
def sources():
    point = Point((0.0,0.0,0.0))
    uniform = spherical_uniform(r_outer=8.85)
    cf252_spec = Watt(1.0250e6,2.926e-6)
    u235_spec = Watt(0.7747e6,4.852e-6)
    u238_spec = Watt(0.6483e6,6.811e-6)
    cf252 = openmc.IndependentSource(space=point,
                                     energy=cf252_spec,
                                     strength=100)
    u235 = openmc.IndependentSource(space=uniform,
                                     energy=u235_spec,
                                     strength=0.5)
    u238 = openmc.IndependentSource(space=uniform,
                                     energy=u238_spec,
                                     strength=57.3)
    return [cf252,u235,u238]
                                                                

@pytest.fixture()
def sphere_model(sources):
    
    mat = openmc.Material(name='heu')
    mat.add_element('U',1,percent_type="ao",enrichment=92.2)
    mat.set_density('g/cm3',18.6)
    
    materials = openmc.Materials([mat])
    
    
    
    surface = openmc.Sphere(r=8.85,boundary_type='vacuum')
    cell = openmc.Cell(fill=mat,region=-surface)
    
    root = openmc.Universe()
    root.add_cell(cell)
    
    geometry = openmc.Geometry(root)
    
    model = openmc.Model(geometry=geometry,materials=materials)
    
    # Settings
    model.settings.particles = 100000
    model.settings.batches = 20
    model.settings.inactive = 10
    model.settings.source = sources
        
    return model


def test_source_effectiveness(sphere_model):
    harness = SourceTestHarness("statepoint.20.h5", model=sphere_model)
    harness.main()
