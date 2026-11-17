"""Alternating Li/Be slabs in a cube; geometry and beam settings are editable."""
import numpy as np
import openmc
from . import config as c


def layer_layout():
    """Return (left, right, material name) for alternating slabs, in cm.

    Lithium is first at the negative face of LAYER_AXIS, then beryllium.
    A final slab is clipped to the cube boundary if the widths do not divide
    the side exactly. Shared surfaces will prevent gaps and overlapping cells.
    """
    sizes=(c.CUBE_SIDE_CM,c.LI_LAYER_THICKNESS_CM,c.BE_LAYER_THICKNESS_CM)
    if not np.all(np.isfinite(sizes)) or min(sizes)<=0:
        raise ValueError('Cube side and both layer widths must be positive and finite')
    if c.LAYER_AXIS not in ('x','y','z'):
        raise ValueError('LAYER_AXIS must be x, y or z')
    half=c.CUBE_SIDE_CM/2
    layers=[];left=-half
    while left<half-1e-12:
        lithium=len(layers)%2==0
        width=c.LI_LAYER_THICKNESS_CM if lithium else c.BE_LAYER_THICKNESS_CM
        right=min(left+width,half)
        layers.append((left,right,'Li' if lithium else 'Be'))
        left=right
    return layers


def source_point():
    """Choose a beam point inside the vacuum world and outside the target.

    For a beam parallel to the slabs, choose the centre of a specified slab;
    the default is a lithium slab near the middle. For a beam normal to the
    slabs, enter the cube at the centre of its upstream face.
    """
    if c.SOURCE_POSITION_CM is not None:
        return tuple(c.SOURCE_POSITION_CM)
    direction=np.asarray(c.SOURCE_DIRECTION,float)
    if direction.shape!=(3,) or not np.all(np.isfinite(direction)):
        raise ValueError('SOURCE_DIRECTION must have three finite coordinates')
    active=np.flatnonzero(direction!=0)
    if len(active)!=1:
        raise ValueError('Automatic source placement needs an axis-aligned direction; otherwise set SOURCE_POSITION_CM')
    if not 0<c.SOURCE_STANDOFF_CM<c.WORLD_MARGIN_CM:
        raise ValueError('Source standoff must be positive and smaller than WORLD_MARGIN_CM')
    beam_axis=int(active[0]);layer_axis='xyz'.index(c.LAYER_AXIS)
    point=np.zeros(3)
    if beam_axis!=layer_axis:
        layers=layer_layout()
        index=c.BEAM_ENTRY_LAYER
        if index is None:
            index=min((j for j,l in enumerate(layers) if l[2]=='Li'),
                      key=lambda j:abs((layers[j][0]+layers[j][1])/2))
        if not isinstance(index,int) or not 0<=index<len(layers):
            raise ValueError('BEAM_ENTRY_LAYER must be a valid zero-based slab index')
        point[layer_axis]=(layers[index][0]+layers[index][1])/2
    point[beam_axis]=-np.sign(direction[beam_axis])*(c.CUBE_SIDE_CM/2+c.SOURCE_STANDOFF_CM)
    return tuple(point)


def build_model(particles=1000,batches=10):
    """Construct a fresh fixed-source cube with an outer vacuum box."""
    layers=layer_layout()
    if not np.isfinite(c.WORLD_MARGIN_CM) or c.WORLD_MARGIN_CM<=0:
        raise ValueError('WORLD_MARGIN_CM must be positive and finite')
    li=openmc.Material(name='Lithium')
    li.set_density('g/cm3',c.LI_DENSITY_G_CM3);li.temperature=c.TEMPERATURE_K
    if c.LI6_ATOM_FRACTION is None:
        li.add_element('Li',1.,percent_type='ao')
    else:
        f=c.LI6_ATOM_FRACTION
        if not 0<=f<=1:raise ValueError('LI6_ATOM_FRACTION must be in [0,1]')
        if f>0:li.add_nuclide('Li6',f)
        if f<1:li.add_nuclide('Li7',1-f)
    be=openmc.Material(name='Beryllium')
    be.add_nuclide('Be9',1.)
    be.set_density('g/cm3',c.BE_DENSITY_G_CM3);be.temperature=c.TEMPERATURE_K
    half=c.CUBE_SIDE_CM/2
    # The target surfaces transmit particles; only the surrounding world is vacuum.
    box=openmc.model.RectangularParallelepiped(-half,half,-half,half,-half,half)
    extent=half+c.WORLD_MARGIN_CM
    world=openmc.model.RectangularParallelepiped(-extent,extent,-extent,extent,-extent,extent,
                                                boundary_type='vacuum')
    plane={'x':openmc.XPlane,'y':openmc.YPlane,'z':openmc.ZPlane}[c.LAYER_AXIS]
    keyword=c.LAYER_AXIS+'0'
    boundaries=[plane(**{keyword:layers[0][0]})]
    boundaries.extend(plane(**{keyword:right}) for _,right,_ in layers)
    cells=[]
    for j,(_,_,name) in enumerate(layers):
        cells.append(openmc.Cell(name=f'Layer {j:02d} {name}',fill=li if name=='Li' else be,
            region=-box & +boundaries[j] & -boundaries[j+1]))
    cells.append(openmc.Cell(name='Surrounding vacuum',region=-world & +box))
    model=openmc.Model(geometry=openmc.Geometry(cells))
    model.settings.run_mode='fixed source'
    model.settings.particles=particles;model.settings.batches=batches
    model.settings.verbosity=1
    model.settings.source=openmc.IndependentSource(
        space=openmc.stats.Point(source_point()),
        angle=openmc.stats.Monodirectional(c.SOURCE_DIRECTION),
        energy=openmc.stats.Discrete([c.SOURCE_ENERGY_EV],[1.]))
    model.settings.output={'summary':False,'tallies':False}
    model.settings.statepoint={'batches':[]}
    return model
