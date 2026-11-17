"""D-shaped TF winding/case envelopes and individual axisymmetric PF coils."""
import numpy as np
import openmc
from ..plasma import contour
from .primitives import revolved
from .central_column import annulus


def make_magnets(c,m,registry):
    base=revolved(contour(c,minor=c.tf_minor_radius,elongation=c.tf_elongation,triangularity=c.tf_triangularity))
    def off(value):return base.offset(float(value))
    outer=off(c.tf_radial_build)
    if outer.points[:,0].min()<=c.column_shield_radii[1] or c.tf_half_width>=outer.points[:,0].min()*np.tan(np.pi/c.tf_coils):
        raise ValueError('TF coils intersect the central shield or exceed their sector width; increase clearances or change coil build')
    if max(np.linalg.norm(outer.points,axis=1))>=c.world_radius:raise ValueError('TF coil exceeds world boundary')
    inset=c.tf_case;ins=c.tf_insulation
    full=(-outer & ~(-base))&+openmc.YPlane(y0=-c.tf_half_width)&-openmc.YPlane(y0=c.tf_half_width)
    insulated=(-off(c.tf_radial_build-inset)&~(-off(inset)))&+openmc.YPlane(y0=-c.tf_half_width+inset)&-openmc.YPlane(y0=c.tf_half_width-inset)
    winding=(-off(c.tf_radial_build-inset-ins)&~(-off(inset+ins)))&+openmc.YPlane(y0=-c.tf_half_width+inset+ins)&-openmc.YPlane(y0=c.tf_half_width-inset-ins)
    # Sector root clipping creates exactly one of the 16 discrete TF coils.
    envelopes=[full]
    cells=[registry.cell('TF winding pack',winding,m['wp'],component='TF_winding'),
           registry.cell('TF insulation',insulated&~winding,m['glass'],component='TF_insulation'),
           registry.cell('TF structural case',full&~insulated,m['ss'],component='TF_case')]
    tf_envelope=full
    for j,p in enumerate(c.pf_coils):
        if np.hypot(p.radius+p.width/2,abs(p.z)+p.height/2)>=c.world_radius:raise ValueError('PF coil exceeds world boundary')
        for r in np.linspace(p.radius-p.width/2,p.radius+p.width/2,4):
            for z in np.linspace(p.z-p.height/2,p.z+p.height/2,4):
                if (r,0,z) in tf_envelope:raise ValueError('PF envelope intersects a TF coil; adjust PF location')
        full=annulus(p.radius-p.width/2,p.radius+p.width/2,p.z-p.height/2,p.z+p.height/2)
        envelopes.append(full)
        winding=annulus(p.radius-p.width/2+c.pf_case,p.radius+p.width/2-c.pf_case,p.z-p.height/2+c.pf_case,p.z+p.height/2-c.pf_case)
        cells.extend([registry.cell(f'PF{j+1} winding',winding,m['wp'],component='PF_winding',pf_id=j+1),
          registry.cell(f'PF{j+1} case',full&~winding,m['ss'],component='PF_case',pf_id=j+1)])
    return cells,envelopes
