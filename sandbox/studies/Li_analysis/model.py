"""Editable lithium sphere, external pencil beam and verified Li response."""
import openmc
from . import config as c


def build_model(particles=1000,batches=10):
    material=openmc.Material(name='Lithium metal')
    material.set_density('g/cm3',c.DENSITY_G_CM3)
    material.temperature=c.TEMPERATURE_K
    if c.LI6_ATOM_FRACTION is None:
        material.add_element('Li',1.,percent_type='ao')
    else:
        fraction=c.LI6_ATOM_FRACTION
        if not 0<=fraction<=1:raise ValueError('Li6 isotope fraction must lie in [0,1]')
        if fraction>0:material.add_nuclide('Li6',fraction)
        if fraction<1:material.add_nuclide('Li7',1-fraction)
    target=openmc.Sphere(r=c.RADIUS_CM)
    world=openmc.Sphere(r=c.RADIUS_CM+c.WORLD_MARGIN_CM,boundary_type='vacuum')
    model=openmc.Model(geometry=openmc.Geometry([
        openmc.Cell(fill=material,region=-target),
        openmc.Cell(region=+target & -world)]))
    model.settings.run_mode='fixed source';model.settings.particles=particles
    model.settings.batches=batches;model.settings.verbosity=1
    model.settings.source=openmc.IndependentSource(
        space=openmc.stats.Point(c.SOURCE_POSITION_CM if c.SOURCE_POSITION_CM is not None
                                 else (0.,0.,-c.RADIUS_CM-c.SOURCE_STANDOFF_CM)),
        angle=openmc.stats.Monodirectional(c.SOURCE_DIRECTION),
        energy=openmc.stats.Discrete([c.SOURCE_ENERGY_EV],[1.]))
    model.settings.output={'summary':False,'tallies':False}
    model.settings.statepoint={'batches':[]}
    return model
