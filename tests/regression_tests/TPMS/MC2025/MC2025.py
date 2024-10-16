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
    mFuel.set_density('g/cm3', 10.045)
    mFuel.add_nuclide("U234", 6.15169E+18)
    mFuel.add_nuclide("U235", 6.89220E+20)
    mFuel.add_nuclide("U236", 3.16265E+18)
    mFuel.add_nuclide("U238", 2.17103E+22)
    mFuel.add_nuclide("C12", 9.13357E+18)
    mFuel.add_nuclide("N14", 1.04072E+19)
    mFuel.add_nuclide("O16", 4.48178E+22)
    mColors[mFuel] = (208, 52, 52)
    materials.append(mFuel)
    return mFuel, materials, mColors

def make_clad(materials:openmc.Materials, mColors:dict) -> tuple[openmc.Material,openmc.Materials,dict]:
    mClad = openmc.Material(1, "mClad", 300.)
    mClad.set_density('g/cm3', 6.56)
    mClad.add_nuclide("O16", 1.19276E-03, percent_type="wo")
    mClad.add_nuclide("O17", 4.82878E-07, percent_type="wo")
    mClad.add_nuclide("O18", 2.75825E-06, percent_type="wo")
    mClad.add_nuclide("Cr50", 4.16117E-05, percent_type="wo")
    mClad.add_nuclide("Cr52", 8.34483E-04, percent_type="wo")
    mClad.add_nuclide("Cr53", 9.64457E-05, percent_type="wo")
    mClad.add_nuclide("Cr54", 2.44600E-05, percent_type="wo")
    mClad.add_nuclide("Fe54", 1.12572E-04, percent_type="wo")
    mClad.add_nuclide("Fe56", 1.83252E-03, percent_type="wo")
    mClad.add_nuclide("Fe57", 4.30778E-05, percent_type="wo")
    mClad.add_nuclide("Fe58", 5.83334E-06, percent_type="wo")
    mClad.add_nuclide("Zr90", 4.97862E-01, percent_type="wo")
    mClad.add_nuclide("Zr91", 1.09780E-01, percent_type="wo")
    mClad.add_nuclide("Zr92", 1.69646E-01, percent_type="wo")
    mClad.add_nuclide("Zr94", 1.75665E-01, percent_type="wo")
    mClad.add_nuclide("Zr96", 2.89038E-02, percent_type="wo")
    mClad.add_nuclide("Sn112", 1.27604E-04, percent_type="wo")
    mClad.add_nuclide("Sn114", 8.83732E-05, percent_type="wo")
    mClad.add_nuclide("Sn115", 4.59255E-05, percent_type="wo")
    mClad.add_nuclide("Sn116", 1.98105E-03, percent_type="wo")
    mClad.add_nuclide("Sn117", 1.05543E-03, percent_type="wo")
    mClad.add_nuclide("Sn118", 3.35688E-03, percent_type="wo")
    mClad.add_nuclide("Sn119", 1.20069E-03, percent_type="wo")
    mClad.add_nuclide("Sn120", 4.59220E-03, percent_type="wo")
    mClad.add_nuclide("Sn122", 6.63497E-04, percent_type="wo")
    mClad.add_nuclide("Sn124", 8.43355E-04, percent_type="wo")
    materials.append(mClad)
    mColors[mClad] = (163, 163, 163)
    return mClad, materials, mColors

def make_cool(materials:openmc.Materials, mColors:dict) -> tuple[openmc.Material,openmc.Materials,dict]:
    mCool = openmc.Material(2, "mCool", 300.)
    mCool.set_density('g/cm3', 0.76973)
    mCool.add_nuclide("H1", 2/3)
    mCool.add_nuclide("O16", 1/3)
    mCool.add_s_alpha_beta("c_H_in_H2O")
    materials.append(mCool)
    mColors[mCool] = (12, 12, 252)
    return mCool, materials, mColors

openmc.Materials.cross_sections = "/Users/fernpa/openmc/libs/endfb80_hdf5/cross_sections.xml"
materials = openmc.Materials()
mColors = {}
mFuel, materials, mColors = make_fuel(materials, mColors)
mClad, materials, mColors = make_clad(materials, mColors)
mCool, materials, mColors = make_cool(materials, mColors)
materials.export_to_xml()

region1 = Box(1.00, "reflective")
tpms1 = openmc.TPMS("Schwarz_P", cst=+0.1, pitch=1.)
tpms2 = openmc.TPMS("Schwarz_P", cst=-0.1, pitch=1.)

cell1 = openmc.Cell(0, "cFuel", mFuel, region1 & +tpms1)
cell2 = openmc.Cell(1, "cClad", mClad, region1 & +tpms2 & -tpms1)
cell3 = openmc.Cell(2, "cCool", mCool, region1 & -tpms2)
universe = openmc.Universe(0, "universe", [cell1,cell2,cell3])
root = openmc.Cell(3, "root", universe, region1)
geometry = openmc.Geometry([root])
geometry.export_to_xml()

settings = openmc.Settings()
settings.batches = 1000
settings.particles = 1000
settings.inactive = 10
settings.seed = 10001
settings.export_to_xml()

plots = openmc.Plots()
plot = openmc.Plot.from_geometry(geometry)
plot.basis = 'xy'
plot.colors = mColors
plot.color_by = 'material'
plot.origin = (0.,0.,0.)
plot.width = (1., 1.)
plot.pixels = (400, 400)
plots.append(plot)
plots.export_to_xml()

# plots = openmc.Plots()
# plot = openmc.Plot.from_geometry(geometry)
# plot.type = 'voxel'
# plot.colors = mColors
# plot.color_by = 'material'
# plot.origin = (0.,0.,0.)
# plot.width = (1., 1., 1.)
# plot.pixels = (400, 400, 400)
# plots.append(plot)
# plots.export_to_xml()