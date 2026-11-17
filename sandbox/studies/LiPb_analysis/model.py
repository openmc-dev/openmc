"""Editable lithium sphere, external pencil beam and verified Li response."""
import openmc
from . import config as c


def build_model(particles=1000,batches=10):
    material = openmc.Material(name="PbLi")
    material.add_element("Li", (1-0.843), percent_type="ao")
    material.add_element("Pb", 0.843, percent_type="ao")
    material.set_density("g/cm3", 9.4)  # Illustrative: set for your operating temperature

    material.temperature=c.TEMPERATURE_K

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
