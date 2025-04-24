import openmc
import numpy as np

import openmc.stats
import openmc.source
import openmc.volume

def Box(c1, c2, c3, boundary_type="transmission"):
    xP0 = openmc.XPlane(-c1/2, boundary_type=boundary_type)
    xP1 = openmc.XPlane(+c1/2, boundary_type=boundary_type)
    yP0 = openmc.YPlane(-c2/2, boundary_type=boundary_type)
    yP1 = openmc.YPlane(+c2/2, boundary_type=boundary_type)
    zP0 = openmc.ZPlane(-c3/2, boundary_type=boundary_type)
    zP1 = openmc.ZPlane(+c3/2, boundary_type=boundary_type)
    return +xP0 & (-xP1 & (+yP0 & (-yP1 & (+zP0 & -zP1))))

def make_fuel(materials:openmc.Materials, mColors:dict) -> tuple[openmc.Material,openmc.Materials,dict]:
    mFuel = openmc.Material(0, "mFuel", 900.)
    mFuel.set_density('g/cm3', 10.5)
    mFuel.add_nuclide("U234", 0.0013)
    mFuel.add_nuclide("U235", 0.1975)
    mFuel.add_nuclide("U236", 0.0010)
    mFuel.add_nuclide("U238", 0.8002)
    mFuel.add_nuclide("O16", 1.9952)
    mFuel.add_nuclide("O17", 0.00074)
    mFuel.add_nuclide("O18", 0.00408)
    mColors[mFuel] = (208, 52, 52)
    materials.append(mFuel)
    return mFuel, materials, mColors

def make_cool(materials:openmc.Materials, mColors:dict) -> tuple[openmc.Material,openmc.Materials,dict]:
    mCool = openmc.Material(2, "mCool", 600.)
    mCool.set_density('g/cm3', 0.110)
    mCool.add_nuclide("C12", 2.69751E-01, "wo")
    mCool.add_nuclide("C13", 3.16149E-03, "wo")
    mCool.add_nuclide("O16", 7.25118E-01, "wo")
    mCool.add_nuclide("O17", 2.93558E-04, "wo")
    mCool.add_nuclide("O18", 1.67683E-03, "wo")
    materials.append(mCool)
    mColors[mCool] = (12, 12, 252)
    return mCool, materials, mColors

def make_refl(materials:openmc.Materials, mColors:dict) -> tuple[openmc.Material,openmc.Materials,dict]:
    mRefl = openmc.Material(3, "mRefl", 600.)
    mRefl.set_density('g/cm3', 0.129)
    mRefl.add_nuclide("C12", 2.69751E-01, "wo")
    mRefl.add_nuclide("C13", 3.16149E-03, "wo")
    mRefl.add_nuclide("O16", 7.25118E-01, "wo")
    mRefl.add_nuclide("O17", 2.93558E-04, "wo")
    mRefl.add_nuclide("O18", 1.67683E-03, "wo")
    materials.append(mRefl)
    mColors[mRefl] = (12, 12, 12)
    return mRefl, materials, mColors

def fPitch(x, y, z):
    ''' Achtung! arguments to pitch function can only be (x,y,z)'''
    # radial pitch radially varying
    pitch = 1-np.sqrt(x**2+y**2)/5
    # pitch = np.full_like(x, 1.)
    return pitch

def isovalue1(x, y, z):
    # thickness = np.full_like(x, +0.1)
    thickness = +0.1 - 0.02*(y+1.5)
    isovalue = 2*np.pi*thickness/fPitch(x,y,z)
    return isovalue

def isovalue2(x, y, z):
    # thickness = np.full_like(x, -0.1)
    thickness = -0.1 + 0.02*(y+1.5)
    isovalue = 2*np.pi*thickness/fPitch(x,y,z)
    return isovalue

materials = openmc.Materials()
mColors = {}
mFuel, materials, mColors = make_fuel(materials, mColors)
mCool, materials, mColors = make_cool(materials, mColors)
mRefl, materials, mColors = make_refl(materials, mColors)
materials.export_to_xml()

box = Box(3., 3., 1., boundary_type="reflective")
tpms1 = openmc.FunctionTPMS.from_interpolated_functions("Schwarz_P", isovalue1, fPitch, (-1.5,+1.5), (-1.5,+1.5), (-0.5,+0.5), (10,10,1))
tpms2 = openmc.FunctionTPMS.from_interpolated_functions("Schwarz_P", isovalue2, fPitch, (-1.5,+1.5), (-1.5,+1.5), (-0.5,+0.5), (10,10,1))


cell1 = openmc.Cell(0, "cFuel", mFuel, box & (-tpms1 & +tpms2))
cell2 = openmc.Cell(2, "cCool", mCool, box & (+tpms1 | -tpms2))

universe = openmc.Universe(0, "universe", [cell1, cell2])

root = openmc.Cell(100, "root", universe, box)
geometry = openmc.Geometry([root])
geometry.export_to_xml()

space = openmc.stats.Point((0.,0.,0.))
source = openmc.IndependentSource(space)
volCalc = openmc.volume.VolumeCalculation([cell1, cell2], 1000000000)
settings = openmc.Settings()
settings.batches = 1000
settings.particles = 10000
settings.inactive = 10
settings.volume_calculations = [volCalc]
settings.source = source
settings.export_to_xml()

plots = openmc.Plots()
plot = openmc.Plot.from_geometry(geometry)
plot.basis = 'xy'
plot.colors = mColors
plot.color_by = 'material'
plot.origin = (0.,0.,0.)
plot.width = (3., 3.)
plot.pixels = (600, 600)
plots.append(plot)
plots.export_to_xml()

plot = openmc.Plot.from_geometry(geometry)
plot.basis = 'xz'
plot.colors = mColors
plot.color_by = 'material'
plot.origin = (0.,0.,0.)
plot.width = (3., 1.)
plot.pixels = (600, 200)
plots.append(plot)
plots.export_to_xml()