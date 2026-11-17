"""Common nested offset contours keep radial interfaces coincident."""
import numpy as np
import openmc
from ..plasma import contour
from .primitives import revolved


def contours(c):
    plasma=revolved(contour(c))
    theta=np.linspace(0,2*np.pi,c.contour_points,endpoint=False)
    build=c.blanket_inboard+(c.blanket_outboard-c.blanket_inboard)*(1+np.cos(theta))/2
    cache={}
    def offset(depth):
        values=np.asarray(depth)+np.zeros(c.contour_points)
        key=tuple(np.round(values,10))
        if key not in cache:cache[key]=plasma.offset(values)
        return cache[key]
    fw=offset(c.plasma_wall_gap)
    blanket_outer=offset(c.plasma_wall_gap+build)
    vv_inner=offset(c.plasma_wall_gap+build+c.vessel_gap)
    vv_wall1=offset(c.plasma_wall_gap+build+c.vessel_gap+c.vessel_inner_wall)
    vv_filler=offset(c.plasma_wall_gap+build+c.vessel_gap+c.vessel_inner_wall+c.vessel_interwall)
    vv_outer=offset(c.plasma_wall_gap+build+c.vessel_gap+c.vessel_inner_wall+c.vessel_interwall+c.vessel_outer_wall)
    shield_inner=offset(c.plasma_wall_gap+build+c.vessel_gap+c.vessel_inner_wall+c.vessel_interwall+c.vessel_outer_wall+c.thermal_gap)
    shield_outer=offset(c.plasma_wall_gap+build+c.vessel_gap+c.vessel_inner_wall+c.vessel_interwall+c.vessel_outer_wall+c.thermal_gap+c.thermal_shield)
    return dict(plasma=plasma,fw=fw,blanket_outer=blanket_outer,vv_inner=vv_inner,
      vv_wall1=vv_wall1,vv_filler=vv_filler,vv_outer=vv_outer,shield_inner=shield_inner,shield_outer=shield_outer,
      offset=offset,build=build)


def make_vessel(c,m,g,registry):
    cells=[]
    for name,inner,outer,material in [('VV inner wall','vv_inner','vv_wall1',m['ss']),
      ('VV interwall shielding coolant','vv_wall1','vv_filler',m['vessel_mix']),
      ('VV outer wall','vv_filler','vv_outer',m['ss']),
      ('Thermal shield','shield_inner','shield_outer',m['ss'])]:
        cells.append(registry.cell(name,-g[outer] & ~(-g[inner]),material,component='vessel'))
    return cells
