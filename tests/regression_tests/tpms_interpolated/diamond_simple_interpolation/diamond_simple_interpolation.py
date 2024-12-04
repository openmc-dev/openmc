import openmc

def Box(c:float, boundary_type="transmission"):
    xP0 = openmc.XPlane(-c/2, boundary_type=boundary_type)
    xP1 = openmc.XPlane(+c/2, boundary_type=boundary_type)
    yP0 = openmc.YPlane(-c/2, boundary_type=boundary_type)
    yP1 = openmc.YPlane(+c/2, boundary_type=boundary_type)
    zP0 = openmc.ZPlane(-c/2, boundary_type=boundary_type)
    zP1 = openmc.ZPlane(+c/2, boundary_type=boundary_type)
    return +xP0 & (-xP1 & (+yP0 & (-yP1 & (+zP0 & -zP1))))

def make_fuel(materials:openmc.Materials, mColors:dict) -> tuple[openmc.Material,openmc.Materials,dict]:
    mFuel = openmc.Material(0, "mFuel", 300.)
    mFuel.set_density('g/cm3', 10.5)
    mFuel.add_nuclide("U235", 0.935/3)
    mFuel.add_nuclide("U238", 0.065/3)
    mFuel.add_nuclide("O16",  2.000/3)
    mColors[mFuel] = (208, 52, 52)
    materials.append(mFuel)
    return mFuel, materials, mColors

def make_cool(materials:openmc.Materials, mColors:dict) -> tuple[openmc.Material,openmc.Materials,dict]:
    mCool = openmc.Material(2, "mCool", 300.)
    mCool.set_density('g/cm3', 0.995)
    mCool.add_nuclide("H1", 2/3)
    mCool.add_nuclide("O16", 1/3)
    mCool.add_s_alpha_beta("c_H_in_H2O")
    materials.append(mCool)
    mColors[mCool] = (12, 12, 252)
    return mCool, materials, mColors

def fPitch(x, y, z):
    r = (x**2+y**2)**0.5
    return 1. + 0.1*r

def fThickness(x, y, z):
    return 0. + 0.05*z**2

materials = openmc.Materials()
mColors = {}
mFuel, materials, mColors = make_fuel(materials, mColors)
mCool, materials, mColors = make_cool(materials, mColors)
materials.export_to_xml()

region1 = Box(5.00, "reflective")
tpms1 = openmc.FunctionTPMS.from_interpolated_functions("Diamond", fThickness, fPitch, (-2.5,+2.5), (-2.5,+2.5), (-2.5,+2.5))

cell1 = openmc.Cell(0, "cFuel", mFuel, region1 & +tpms1) # good practice : always put the tpms in the box
cell2 = openmc.Cell(2, "cCool", mCool, region1 & -tpms1)
universe = openmc.Universe(0, "universe", [cell1,cell2])
root = openmc.Cell(3, "root", universe, region1)
geometry = openmc.Geometry([root])
geometry.export_to_xml()

settings = openmc.Settings()
settings.batches = 100
settings.particles = 100
settings.inactive = 10
settings.export_to_xml()

plots = openmc.Plots()
plot = openmc.Plot.from_geometry(geometry)
plot.basis = 'xy'
plot.colors = mColors
plot.color_by = 'material'
plot.origin = (0.,0.,0.)
plot.width = (5., 5.)
plot.pixels = (600, 600)
plots.append(plot)
plots.export_to_xml()

plot = openmc.Plot.from_geometry(geometry)
plot.basis = 'xz'
plot.colors = mColors
plot.color_by = 'material'
plot.origin = (0.,0.,0.)
plot.width = (5., 5.)
plot.pixels = (600, 600)
plots.append(plot)
plots.export_to_xml()

plot = openmc.Plot.from_geometry(geometry)
plot.type = 'voxel'
plot.colors = mColors
plot.color_by = 'material'
plot.origin = (0.,0.,0.)
plot.width = (5., 5., 5.)
plot.pixels = (600, 600, 600)
plots.append(plot)
plots.export_to_xml()