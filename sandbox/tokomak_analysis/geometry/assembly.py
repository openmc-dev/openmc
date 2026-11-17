"""Full 360-degree CSG assembled as 16 repeated native universes."""
from dataclasses import dataclass
from time import perf_counter
import numpy as np
import openmc
from ..config import Config
from ..materials import make_materials
from .primitives import Registry,sector,union,share_interfaces
from .vessel import contours,make_vessel
from .central_column import make_column
from .magnets import make_magnets
from .divertor import make_divertor
from .blanket import make_blanket


@dataclass
class Assembly:
    model: openmc.Model
    config: Config
    registry: Registry
    contours: dict
    layers: list
    construction_seconds: float


def build_model(config=None,*,tallies=True,output=None):
    c=(config or Config()).validate();start=perf_counter();openmc.reset_auto_ids()
    m=make_materials(c);g=contours(c);registry=Registry()
    world=openmc.Sphere(r=c.world_radius,boundary_type='vacuum')
    plasma=registry.cell('Plasma source domain',-g['plasma'],component='plasma')
    blanket,layers=make_blanket(c,m,g,registry)
    divertor,divertor_envelopes=make_divertor(c,m,registry)
    vessel=make_vessel(c,m,g,registry)
    column,column_envelopes=make_column(c,m,registry)
    magnets,magnet_envelopes=make_magnets(c,m,registry)
    inside=[plasma,blanket,*divertor]
    interior_vacuum=registry.cell('SOL chamber and connected divertor vacuum',-g['vv_inner']&~union([plasma.region,blanket.region,*divertor_envelopes]),component='chamber_vacuum')
    occupied=[*inside,interior_vacuum,*vessel,*column,*magnets]
    sector_universe=openmc.Universe(cells=occupied,name='DEMO repeated sector')
    external=registry.cell('External cryostat vacuum',-world&~union([-g['vv_outer'],vessel[-1].region,*column_envelopes,*magnet_envelopes]),component='external_vacuum')
    sector_universe.add_cell(external)
    root=openmc.Universe();width=2*np.pi/c.tf_coils
    for j in range(c.tf_coils):
        centre=j*width
        parent=registry.cell(f'Reactor sector {j}',-world&sector(centre-width/2,centre+width/2),sector_universe,component='sector',sector=j)
        parent.rotation=(0.,0.,np.rad2deg(centre));root.add_cell(parent)
    geometry=openmc.Geometry(root);share_interfaces(geometry)
    model=openmc.Model(geometry=geometry)
    model.settings.run_mode='fixed source';model.settings.particles=c.particles;model.settings.batches=c.batches;model.settings.seed=c.seed
    model.settings.temperature={'method':'nearest','tolerance':c.temperature_tolerance}
    model.settings.photon_transport=c.coupled_photons
    model.settings.output={'summary':True,'tallies':False}
    if output is not None:
        from pathlib import Path
        model.settings.output['path']=str(Path(output).resolve())+'/'
    model.settings.statepoint={'batches':[c.batches]}
    model.settings.verbosity=5
    model.settings.max_lost_particles=1
    model.settings.rel_max_lost_particles=1e-15
    from ..neutronics.source import make_source
    model.settings.source=make_source(c,plasma,g['plasma'].points,output)
    assembly=Assembly(model,c,registry,g,layers,perf_counter()-start)
    if tallies:
        from ..neutronics.tallies import add_tallies
        add_tallies(assembly)
    return assembly
