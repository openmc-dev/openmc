"""Lower single-null-like cassette and inclined target stack approximations."""
import numpy as np
from .primitives import revolved,rz_halfspace


def target_strip(a,b,lo,hi,sign):
    a=np.asarray(a);b=np.asarray(b);t=(b-a)/np.linalg.norm(b-a)
    normal=sign*np.array([-t[1],t[0]])
    return (rz_halfspace(t,np.dot(t,a)) & rz_halfspace(-t,-np.dot(t,b))
            & rz_halfspace(normal,np.dot(normal,a)+lo)
            & rz_halfspace(-normal,-np.dot(normal,a)-hi))


def rectangle(rlo,rhi,zlo,zhi):return -revolved([(rlo,zlo),(rhi,zlo),(rhi,zhi),(rlo,zhi)])


def make_divertor(c,m,registry):
    # Coordinates tied to plasma major radius/height. Geometry checks enforce
    # compatibility with the VV and plasma for a changed radial build.
    def p(r,z):return (c.major_radius+r,-c.elongation*c.minor_radius+z)
    inner=[p(*q) for q in c.divertor_inner_target];outer_target=[p(*q) for q in c.divertor_outer_target]
    cells=[];envelopes=[]
    for label,a,b,sign in [('Inner',*inner,-1),('Outer',*outer_target,1)]:
        previous=None
        lo=0.
        for layer,width,key in [('W armour',c.divertor_armour,'w'),('CuCrZr heat sink',c.divertor_copper,'cucrzr'),('water cooling',c.divertor_water,'water'),('steel support',c.divertor_support,'steel')]:
            hi=lo+width
            cumulative=target_strip(a,b,0.,hi,sign)
            region=cumulative if previous is None else cumulative&~previous
            cells.append(registry.cell(label+' divertor '+layer,region,m[key],component='divertor'))
            previous=cumulative;lo=hi
        envelopes.append(previous)
    r0,r1=[c.major_radius+x for x in c.cassette_radial_offsets]
    z0,z1=[-c.elongation*c.minor_radius+x for x in c.cassette_vertical_offsets]
    cr0,cr1=[c.major_radius+x for x in c.cassette_coolant_radial_offsets]
    cz0,cz1=[-c.elongation*c.minor_radius+x for x in c.cassette_coolant_vertical_offsets]
    outer=rectangle(r0,r1,z0,z1)
    coolant=rectangle(cr0,cr1,cz0,cz1)
    cells.extend([registry.cell('Divertor cassette',outer&~coolant,m['steel'],component='divertor'),
      registry.cell('Divertor cassette cooling',coolant,m['water'],component='divertor'),
      registry.cell('Divertor shielding liner',rectangle(r0,r1,z1,z1+c.cassette_liner),m['w'],component='divertor')])
    return cells,[*envelopes,outer,rectangle(r0,r1,z1,z1+c.cassette_liner)]
