import openmc
import pytest
import os

# compile the external source
def compile_source():
    # needs a more robust way to know where the
    # Makefile is
    status  = os.system('cp /usr/local/share/Makefile .')
    assert os.WEXITSTATUS(status) == 0
    status = os.system('make')
    assert os.WEXITSTATUS(status) == 0
    
    return

# build the test geometry
def make_geometry():
    mats = openmc.Materials()

    natural_lead = openmc.Material(1, "natural_lead")
    natural_lead.add_element('Pb', 1,'ao')
    natural_lead.set_density('g/cm3', 11.34)
    mats.append(natural_lead)

    # surfaces
    surface_sph1 = openmc.Sphere(r=100) 
    volume_sph1 = -surface_sph1 
    
    # cell
    cell_1 = openmc.Cell(region=volume_sph1)
    cell_1.fill = natural_lead #assigning a material to a cell
    universe = openmc.Universe(cells=[cell_1]) #hint, this list will need to include the new cell
    geom = openmc.Geometry(universe)

    # settings
    sett = openmc.Settings()
    batches = 10
    sett.batches = batches
    sett.inactive = 0
    sett.particles = 1000
    sett.particle = "neutron"
    sett.run_mode = 'fixed source'

    #source
    source = openmc.Source()
    source.library = PATH+'source_sampling.so'
    sett.source = source

    # run
    model = openmc.model.Model(geom,mats,sett)
    return model

def test_dlopen_source():
    compile_source()
    model = make_geometry()
    model.run()
    
