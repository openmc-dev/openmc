"""
ICT014-1: RA-6
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Al-U3Si2
mat1 = openmc.Material(material_id=1, name="Al-U3Si2")
mat1.set_density("atom/b-cm", 5.091724e-02)
mat1.add_nuclide("U234", 1.7680e-5)
mat1.add_nuclide("U235", 2.4170e-3)
mat1.add_nuclide("U236", 1.3635e-5)
mat1.add_nuclide("U238", 9.6538e-3)
mat1.add_element("Si", 8.3622e-3)
mat1.add_element("Al", 3.0451e-2)
mat1.add_element("B", 1.9296e-6)

# Al6061
mat2 = openmc.Material(material_id=2, name="Al6061")
mat2.set_density("atom/b-cm", 6.004096e-02)
mat2.add_element("Al", 5.8811e-2)
mat2.add_element("Cu", 6.6527e-5)
mat2.add_element("Cr", 3.4398e-5)
mat2.add_element("Mg", 6.6229e-4)
mat2.add_element("Si", 3.8789e-4)
mat2.add_element("Zn", 2.4866e-7)
mat2.add_element("Fe", 7.8610e-5)

# H2O
mat3 = openmc.Material(material_id=3, name="H2O")
mat3.set_density("atom/b-cm", 1.000130e-01)
mat3.add_element("H", 6.6675e-2)
mat3.add_nuclide("O16", 3.3338e-2)
mat3.add_s_alpha_beta("c_H_in_H2O")

# Cd
mat4 = openmc.Material(material_id=4, name="Cd")
mat4.set_density("atom/b-cm", 4.634000e-02)
mat4.add_element("Cd", 4.6340e-2)

# Al
mat5 = openmc.Material(material_id=5, name="Al")
mat5.set_density("atom/b-cm", 6.026200e-02)
mat5.add_element("Al", 6.0262e-2)

# Ag-In-Cd
mat6 = openmc.Material(material_id=6, name="Ag-In-Cd")
mat6.set_density("atom/b-cm", 5.617460e-02)
mat6.add_element("Ag", 4.5365e-2)
mat6.add_element("In", 7.9765e-3)
mat6.add_element("Cd", 2.8331e-3)

# Stainless steel
mat7 = openmc.Material(material_id=7, name="Stainless steel")
mat7.set_density("atom/b-cm", 8.728029e-02)
mat7.add_element("Fe", 5.9899e-2)
mat7.add_element("Ni", 8.1984e-3)
mat7.add_element("Cr", 1.7582e-2)
mat7.add_element("Si", 6.4246e-4)
mat7.add_element("Mn", 8.7583e-4)
mat7.add_element("C", 6.0091e-5)
mat7.add_element("S", 2.2505e-5)
mat7.add_s_alpha_beta("c_Graphite")

# Al2O3
mat8 = openmc.Material(material_id=8, name="Al2O3")
mat8.set_density("atom/b-cm", 1.169440e-01)
mat8.add_element("Al", 4.6778e-2)
mat8.add_nuclide("O16", 7.0166e-2)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# Al6061 side plates: inner
# Rectangular parallelepiped (6 planes)
surf11_xmin = openmc.XPlane(surface_id=10000, x0=-3.3)
surf11_xmax = openmc.XPlane(surface_id=10001, x0=3.3)
surf11_ymin = openmc.YPlane(surface_id=10002, y0=-4)
surf11_ymax = openmc.YPlane(surface_id=10003, y0=4)
surf11_zmin = openmc.ZPlane(surface_id=10004, z0=-38.35)
surf11_zmax = openmc.ZPlane(surface_id=10005, z0=38.35)

# Al6061 side plates: outer
# Rectangular parallelepiped (6 planes)
surf12_xmin = openmc.XPlane(surface_id=10006, x0=-3.8)
surf12_xmax = openmc.XPlane(surface_id=10007, x0=3.8)
surf12_ymin = openmc.YPlane(surface_id=10008, y0=-4)
surf12_ymax = openmc.YPlane(surface_id=10009, y0=4)
surf12_zmin = openmc.ZPlane(surface_id=10010, z0=-38.35)
surf12_zmax = openmc.ZPlane(surface_id=10011, z0=38.35)

# Al-U3Si2 fuel
# Rectangular parallelepiped (6 planes)
surf21_xmin = openmc.XPlane(surface_id=10012, x0=-3)
surf21_xmax = openmc.XPlane(surface_id=10013, x0=3)
surf21_ymin = openmc.YPlane(surface_id=10014, y0=-0.02533)
surf21_ymax = openmc.YPlane(surface_id=10015, y0=0.02533)
surf21_zmin = openmc.ZPlane(surface_id=10016, z0=-31.15)
surf21_zmax = openmc.ZPlane(surface_id=10017, z0=31.15)

# Al6061 clad: inner fuel plates
# Rectangular parallelepiped (6 planes)
surf22_xmin = openmc.XPlane(surface_id=10018, x0=-3.525)
surf22_xmax = openmc.XPlane(surface_id=10019, x0=3.525)
surf22_ymin = openmc.YPlane(surface_id=10020, y0=-0.0745)
surf22_ymax = openmc.YPlane(surface_id=10021, y0=0.0745)
surf22_zmin = openmc.ZPlane(surface_id=10022, z0=-33.55)
surf22_zmax = openmc.ZPlane(surface_id=10023, z0=33.55)

# Water slot:  inner fuel plates
# Rectangular parallelepiped (6 planes)
surf23_xmin = openmc.XPlane(surface_id=10024, x0=-3.55)
surf23_xmax = openmc.XPlane(surface_id=10025, x0=3.55)
surf23_ymin = openmc.YPlane(surface_id=10026, y0=-0.08)
surf23_ymax = openmc.YPlane(surface_id=10027, y0=0.08)
surf23_zmin = openmc.ZPlane(surface_id=10028, z0=-33.85)
surf23_zmax = openmc.ZPlane(surface_id=10029, z0=35.15)

# Al6061 clad: outer fuel plates
# Rectangular parallelepiped (6 planes)
surf24_xmin = openmc.XPlane(surface_id=10030, x0=-3.525)
surf24_xmax = openmc.XPlane(surface_id=10031, x0=3.525)
surf24_ymin = openmc.YPlane(surface_id=10032, y0=-0.0745)
surf24_ymax = openmc.YPlane(surface_id=10033, y0=0.0745)
surf24_zmin = openmc.ZPlane(surface_id=10034, z0=-38.35)
surf24_zmax = openmc.ZPlane(surface_id=10035, z0=35.15)

# Water slot:  outer fuel plates
# Rectangular parallelepiped (6 planes)
surf25_xmin = openmc.XPlane(surface_id=10036, x0=-3.55)
surf25_xmax = openmc.XPlane(surface_id=10037, x0=3.55)
surf25_ymin = openmc.YPlane(surface_id=10038, y0=-0.08)
surf25_ymax = openmc.YPlane(surface_id=10039, y0=0.08)
surf25_zmin = openmc.ZPlane(surface_id=10040, z0=-38.35)
surf25_zmax = openmc.ZPlane(surface_id=10041, z0=35.15)

# Cd wire
surf26 = openmc.ZCylinder(surface_id=26, r=0.02425)

# Al6061 crossbar
surf29 = openmc.XCylinder(surface_id=29, r=0.625)

# Al6061 nozzle: inner
surf30 = openmc.ZCylinder(surface_id=30, r=2.4895)

# Al6061 nozzle: outer
surf31 = openmc.ZCylinder(surface_id=31, r=3.0895)

# Al6061 side plates: inner
# Rectangular parallelepiped (6 planes)
surf41_xmin = openmc.XPlane(surface_id=10042, x0=-3.3)
surf41_xmax = openmc.XPlane(surface_id=10043, x0=3.3)
surf41_ymin = openmc.YPlane(surface_id=10044, y0=-4)
surf41_ymax = openmc.YPlane(surface_id=10045, y0=4)
surf41_zmin = openmc.ZPlane(surface_id=10046, z0=-38.35)
surf41_zmax = openmc.ZPlane(surface_id=10047, z0=39.65)

# Al6061 side plates: outer
# Rectangular parallelepiped (6 planes)
surf42_xmin = openmc.XPlane(surface_id=10048, x0=-3.8)
surf42_xmax = openmc.XPlane(surface_id=10049, x0=3.8)
surf42_ymin = openmc.YPlane(surface_id=10050, y0=-4)
surf42_ymax = openmc.YPlane(surface_id=10051, y0=4)
surf42_zmin = openmc.ZPlane(surface_id=10052, z0=-38.35)
surf42_zmax = openmc.ZPlane(surface_id=10053, z0=39.65)

# Al6061 internal guide plate
# Rectangular parallelepiped (6 planes)
surf43_xmin = openmc.XPlane(surface_id=10054, x0=-3.525)
surf43_xmax = openmc.XPlane(surface_id=10055, x0=3.525)
surf43_ymin = openmc.YPlane(surface_id=10056, y0=-0.065)
surf43_ymax = openmc.YPlane(surface_id=10057, y0=0.065)
surf43_zmin = openmc.ZPlane(surface_id=10058, z0=-33.55)
surf43_zmax = openmc.ZPlane(surface_id=10059, z0=33.55)

# Al6061 external guide plate
# Rectangular parallelepiped (6 planes)
surf44_xmin = openmc.XPlane(surface_id=10060, x0=-3.525)
surf44_xmax = openmc.XPlane(surface_id=10061, x0=3.525)
surf44_ymin = openmc.YPlane(surface_id=10062, y0=-0.065)
surf44_ymax = openmc.YPlane(surface_id=10063, y0=0.065)
surf44_zmin = openmc.ZPlane(surface_id=10064, z0=-38.35)
surf44_zmax = openmc.ZPlane(surface_id=10065, z0=39.65)

# In-Ag-Cd
# Rectangular parallelepiped (6 planes)
surf50_xmin = openmc.XPlane(surface_id=10066, x0=-3.09)
surf50_xmax = openmc.XPlane(surface_id=10067, x0=3.09)
surf50_ymin = openmc.YPlane(surface_id=10068, y0=-0.11)
surf50_ymax = openmc.YPlane(surface_id=10069, y0=0.11)
surf50_zmin = openmc.ZPlane(surface_id=10070, z0=31.55)
surf50_zmax = openmc.ZPlane(surface_id=10071, z0=94.95)

# edge
surf51 = openmc.XPlane(surface_id=51, x0=-2.98)

# edge
surf52 = openmc.XPlane(surface_id=52, x0=2.98)

# edge
surf53 = openmc.ZCylinder(surface_id=53, r=0.11)

# edge
surf54 = openmc.ZCylinder(surface_id=54, r=0.11)

# indentation
# COG surface type "s" with parameters: 0.1625 tr 0 0.1725 63.25
# This surface type requires manual translation to OpenMC
surf55 = openmc.Sphere(surface_id=55, r=1.0)  # PLACEHOLDER - REPLACE THIS

# SS304L cladding: inner
# Rectangular parallelepiped (6 planes)
surf56_xmin = openmc.XPlane(surface_id=10072, x0=-3.175)
surf56_xmax = openmc.XPlane(surface_id=10073, x0=3.175)
surf56_ymin = openmc.YPlane(surface_id=10074, y0=-0.145)
surf56_ymax = openmc.YPlane(surface_id=10075, y0=0.145)
surf56_zmin = openmc.ZPlane(surface_id=10076, z0=31.55)
surf56_zmax = openmc.ZPlane(surface_id=10077, z0=94.95)

# SS304L cladding: outer
# Rectangular parallelepiped (6 planes)
surf57_xmin = openmc.XPlane(surface_id=10078, x0=-3.245)
surf57_xmax = openmc.XPlane(surface_id=10079, x0=3.245)
surf57_ymin = openmc.YPlane(surface_id=10080, y0=-0.215)
surf57_ymax = openmc.YPlane(surface_id=10081, y0=0.215)
surf57_zmin = openmc.ZPlane(surface_id=10082, z0=30.45)
surf57_zmax = openmc.ZPlane(surface_id=10083, z0=117.15)

# BNCT
# Rectangular parallelepiped (6 planes)
surf60_xmin = openmc.XPlane(surface_id=10084, x0=-38.55)
surf60_xmax = openmc.XPlane(surface_id=10085, x0=38.55)
surf60_ymin = openmc.YPlane(surface_id=10086, y0=-199.0)
surf60_ymax = openmc.YPlane(surface_id=10087, y0=-32.40)
surf60_zmin = openmc.ZPlane(surface_id=10088, z0=-40.75)
surf60_zmax = openmc.ZPlane(surface_id=10089, z0=41.60)

# Al/Cd
surf61 = openmc.YPlane(surface_id=61, y0=-49.50)

# Cd/Al
surf62 = openmc.YPlane(surface_id=62, y0=-49.55)

# Al/Cd
surf63 = openmc.YPlane(surface_id=63, y0=-59.55)

# Cd/Al2O3
surf64 = openmc.YPlane(surface_id=64, y0=-59.70)

# Pure aluminum grid plate
# Rectangular parallelepiped (6 planes)
surf90_xmin = openmc.XPlane(surface_id=10090, x0=-30.8)
surf90_xmax = openmc.XPlane(surface_id=10091, x0=30.8)
surf90_ymin = openmc.YPlane(surface_id=10092, y0=-40.5)
surf90_ymax = openmc.YPlane(surface_id=10093, y0=40.5)
surf90_zmin = openmc.ZPlane(surface_id=10094, z0=-60.75)
surf90_zmax = openmc.ZPlane(surface_id=10095, z0=-40.75)

# Primary water hole
surf91 = openmc.ZCylinder(surface_id=91, r=3.0895)

# 1st (outer) plate and slot
# Rectangular parallelepiped (6 planes)
surf101_xmin = openmc.XPlane(surface_id=10096, x0=-3.55)
surf101_xmax = openmc.XPlane(surface_id=10097, x0=3.55)
surf101_ymin = openmc.YPlane(surface_id=10098, y0=-0.08)
surf101_ymax = openmc.YPlane(surface_id=10099, y0=0.08)
surf101_zmin = openmc.ZPlane(surface_id=10100, z0=-38.35)
surf101_zmax = openmc.ZPlane(surface_id=10101, z0=35.15)

# 2nd (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf102_xmin = openmc.XPlane(surface_id=10102, x0=-3.55)
surf102_xmax = openmc.XPlane(surface_id=10103, x0=3.55)
surf102_ymin = openmc.YPlane(surface_id=10104, y0=-0.08)
surf102_ymax = openmc.YPlane(surface_id=10105, y0=0.08)
surf102_zmin = openmc.ZPlane(surface_id=10106, z0=-33.85)
surf102_zmax = openmc.ZPlane(surface_id=10107, z0=35.15)

# 3rd (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf103_xmin = openmc.XPlane(surface_id=10108, x0=-3.55)
surf103_xmax = openmc.XPlane(surface_id=10109, x0=3.55)
surf103_ymin = openmc.YPlane(surface_id=10110, y0=-0.08)
surf103_ymax = openmc.YPlane(surface_id=10111, y0=0.08)
surf103_zmin = openmc.ZPlane(surface_id=10112, z0=-33.85)
surf103_zmax = openmc.ZPlane(surface_id=10113, z0=35.15)

# 4th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf104_xmin = openmc.XPlane(surface_id=10114, x0=-3.55)
surf104_xmax = openmc.XPlane(surface_id=10115, x0=3.55)
surf104_ymin = openmc.YPlane(surface_id=10116, y0=-0.08)
surf104_ymax = openmc.YPlane(surface_id=10117, y0=0.08)
surf104_zmin = openmc.ZPlane(surface_id=10118, z0=-33.85)
surf104_zmax = openmc.ZPlane(surface_id=10119, z0=35.15)

# 5th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf105_xmin = openmc.XPlane(surface_id=10120, x0=-3.55)
surf105_xmax = openmc.XPlane(surface_id=10121, x0=3.55)
surf105_ymin = openmc.YPlane(surface_id=10122, y0=-0.08)
surf105_ymax = openmc.YPlane(surface_id=10123, y0=0.08)
surf105_zmin = openmc.ZPlane(surface_id=10124, z0=-33.85)
surf105_zmax = openmc.ZPlane(surface_id=10125, z0=35.15)

# 6th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf106_xmin = openmc.XPlane(surface_id=10126, x0=-3.55)
surf106_xmax = openmc.XPlane(surface_id=10127, x0=3.55)
surf106_ymin = openmc.YPlane(surface_id=10128, y0=-0.08)
surf106_ymax = openmc.YPlane(surface_id=10129, y0=0.08)
surf106_zmin = openmc.ZPlane(surface_id=10130, z0=-33.85)
surf106_zmax = openmc.ZPlane(surface_id=10131, z0=35.15)

# 7th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf107_xmin = openmc.XPlane(surface_id=10132, x0=-3.55)
surf107_xmax = openmc.XPlane(surface_id=10133, x0=3.55)
surf107_ymin = openmc.YPlane(surface_id=10134, y0=-0.08)
surf107_ymax = openmc.YPlane(surface_id=10135, y0=0.08)
surf107_zmin = openmc.ZPlane(surface_id=10136, z0=-33.85)
surf107_zmax = openmc.ZPlane(surface_id=10137, z0=35.15)

# 8th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf108_xmin = openmc.XPlane(surface_id=10138, x0=-3.55)
surf108_xmax = openmc.XPlane(surface_id=10139, x0=3.55)
surf108_ymin = openmc.YPlane(surface_id=10140, y0=-0.08)
surf108_ymax = openmc.YPlane(surface_id=10141, y0=0.08)
surf108_zmin = openmc.ZPlane(surface_id=10142, z0=-33.85)
surf108_zmax = openmc.ZPlane(surface_id=10143, z0=35.15)

# 9th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf109_xmin = openmc.XPlane(surface_id=10144, x0=-3.55)
surf109_xmax = openmc.XPlane(surface_id=10145, x0=3.55)
surf109_ymin = openmc.YPlane(surface_id=10146, y0=-0.08)
surf109_ymax = openmc.YPlane(surface_id=10147, y0=0.08)
surf109_zmin = openmc.ZPlane(surface_id=10148, z0=-33.85)
surf109_zmax = openmc.ZPlane(surface_id=10149, z0=35.15)

# 11th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf111_xmin = openmc.XPlane(surface_id=10150, x0=-3.55)
surf111_xmax = openmc.XPlane(surface_id=10151, x0=3.55)
surf111_ymin = openmc.YPlane(surface_id=10152, y0=-0.08)
surf111_ymax = openmc.YPlane(surface_id=10153, y0=0.08)
surf111_zmin = openmc.ZPlane(surface_id=10154, z0=-33.85)
surf111_zmax = openmc.ZPlane(surface_id=10155, z0=35.15)

# 12th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf112_xmin = openmc.XPlane(surface_id=10156, x0=-3.55)
surf112_xmax = openmc.XPlane(surface_id=10157, x0=3.55)
surf112_ymin = openmc.YPlane(surface_id=10158, y0=-0.08)
surf112_ymax = openmc.YPlane(surface_id=10159, y0=0.08)
surf112_zmin = openmc.ZPlane(surface_id=10160, z0=-33.85)
surf112_zmax = openmc.ZPlane(surface_id=10161, z0=35.15)

# 13th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf113_xmin = openmc.XPlane(surface_id=10162, x0=-3.55)
surf113_xmax = openmc.XPlane(surface_id=10163, x0=3.55)
surf113_ymin = openmc.YPlane(surface_id=10164, y0=-0.08)
surf113_ymax = openmc.YPlane(surface_id=10165, y0=0.08)
surf113_zmin = openmc.ZPlane(surface_id=10166, z0=-33.85)
surf113_zmax = openmc.ZPlane(surface_id=10167, z0=35.15)

# 14th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf114_xmin = openmc.XPlane(surface_id=10168, x0=-3.55)
surf114_xmax = openmc.XPlane(surface_id=10169, x0=3.55)
surf114_ymin = openmc.YPlane(surface_id=10170, y0=-0.08)
surf114_ymax = openmc.YPlane(surface_id=10171, y0=0.08)
surf114_zmin = openmc.ZPlane(surface_id=10172, z0=-33.85)
surf114_zmax = openmc.ZPlane(surface_id=10173, z0=35.15)

# 15th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf115_xmin = openmc.XPlane(surface_id=10174, x0=-3.55)
surf115_xmax = openmc.XPlane(surface_id=10175, x0=3.55)
surf115_ymin = openmc.YPlane(surface_id=10176, y0=-0.08)
surf115_ymax = openmc.YPlane(surface_id=10177, y0=0.08)
surf115_zmin = openmc.ZPlane(surface_id=10178, z0=-33.85)
surf115_zmax = openmc.ZPlane(surface_id=10179, z0=35.15)

# 16th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf116_xmin = openmc.XPlane(surface_id=10180, x0=-3.55)
surf116_xmax = openmc.XPlane(surface_id=10181, x0=3.55)
surf116_ymin = openmc.YPlane(surface_id=10182, y0=-0.08)
surf116_ymax = openmc.YPlane(surface_id=10183, y0=0.08)
surf116_zmin = openmc.ZPlane(surface_id=10184, z0=-33.85)
surf116_zmax = openmc.ZPlane(surface_id=10185, z0=35.15)

# 17th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf117_xmin = openmc.XPlane(surface_id=10186, x0=-3.55)
surf117_xmax = openmc.XPlane(surface_id=10187, x0=3.55)
surf117_ymin = openmc.YPlane(surface_id=10188, y0=-0.08)
surf117_ymax = openmc.YPlane(surface_id=10189, y0=0.08)
surf117_zmin = openmc.ZPlane(surface_id=10190, z0=-33.85)
surf117_zmax = openmc.ZPlane(surface_id=10191, z0=35.15)

# 18th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf118_xmin = openmc.XPlane(surface_id=10192, x0=-3.55)
surf118_xmax = openmc.XPlane(surface_id=10193, x0=3.55)
surf118_ymin = openmc.YPlane(surface_id=10194, y0=-0.08)
surf118_ymax = openmc.YPlane(surface_id=10195, y0=0.08)
surf118_zmin = openmc.ZPlane(surface_id=10196, z0=-33.85)
surf118_zmax = openmc.ZPlane(surface_id=10197, z0=35.15)

# 19th (outer) plate and slot
# Rectangular parallelepiped (6 planes)
surf119_xmin = openmc.XPlane(surface_id=10198, x0=-3.55)
surf119_xmax = openmc.XPlane(surface_id=10199, x0=3.55)
surf119_ymin = openmc.YPlane(surface_id=10200, y0=-0.08)
surf119_ymax = openmc.YPlane(surface_id=10201, y0=0.08)
surf119_zmin = openmc.ZPlane(surface_id=10202, z0=-38.35)
surf119_zmax = openmc.ZPlane(surface_id=10203, z0=35.15)

# Slot for Cd wire: outer fuel plate
# Rectangular parallelepiped (6 planes)
surf121_xmin = openmc.XPlane(surface_id=10204, x0=-0.03)
surf121_xmax = openmc.XPlane(surface_id=10205, x0=0.03)
surf121_ymin = openmc.YPlane(surface_id=10206, y0=-0.025)
surf121_ymax = openmc.YPlane(surface_id=10207, y0=0.025)
surf121_zmin = openmc.ZPlane(surface_id=10208, z0=-38.35)
surf121_zmax = openmc.ZPlane(surface_id=10209, z0=31.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf123_xmin = openmc.XPlane(surface_id=10210, x0=-0.03)
surf123_xmax = openmc.XPlane(surface_id=10211, x0=0.03)
surf123_ymin = openmc.YPlane(surface_id=10212, y0=-0.025)
surf123_ymax = openmc.YPlane(surface_id=10213, y0=0.025)
surf123_zmin = openmc.ZPlane(surface_id=10214, z0=-33.85)
surf123_zmax = openmc.ZPlane(surface_id=10215, z0=31.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf125_xmin = openmc.XPlane(surface_id=10216, x0=-0.03)
surf125_xmax = openmc.XPlane(surface_id=10217, x0=0.03)
surf125_ymin = openmc.YPlane(surface_id=10218, y0=-0.025)
surf125_ymax = openmc.YPlane(surface_id=10219, y0=0.025)
surf125_zmin = openmc.ZPlane(surface_id=10220, z0=-33.85)
surf125_zmax = openmc.ZPlane(surface_id=10221, z0=31.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf127_xmin = openmc.XPlane(surface_id=10222, x0=-0.03)
surf127_xmax = openmc.XPlane(surface_id=10223, x0=0.03)
surf127_ymin = openmc.YPlane(surface_id=10224, y0=-0.025)
surf127_ymax = openmc.YPlane(surface_id=10225, y0=0.025)
surf127_zmin = openmc.ZPlane(surface_id=10226, z0=-33.85)
surf127_zmax = openmc.ZPlane(surface_id=10227, z0=31.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf129_xmin = openmc.XPlane(surface_id=10228, x0=-0.03)
surf129_xmax = openmc.XPlane(surface_id=10229, x0=0.03)
surf129_ymin = openmc.YPlane(surface_id=10230, y0=-0.025)
surf129_ymax = openmc.YPlane(surface_id=10231, y0=0.025)
surf129_zmin = openmc.ZPlane(surface_id=10232, z0=-33.85)
surf129_zmax = openmc.ZPlane(surface_id=10233, z0=31.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf131_xmin = openmc.XPlane(surface_id=10234, x0=-0.03)
surf131_xmax = openmc.XPlane(surface_id=10235, x0=0.03)
surf131_ymin = openmc.YPlane(surface_id=10236, y0=-0.025)
surf131_ymax = openmc.YPlane(surface_id=10237, y0=0.025)
surf131_zmin = openmc.ZPlane(surface_id=10238, z0=-33.85)
surf131_zmax = openmc.ZPlane(surface_id=10239, z0=31.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf133_xmin = openmc.XPlane(surface_id=10240, x0=-0.03)
surf133_xmax = openmc.XPlane(surface_id=10241, x0=0.03)
surf133_ymin = openmc.YPlane(surface_id=10242, y0=-0.025)
surf133_ymax = openmc.YPlane(surface_id=10243, y0=0.025)
surf133_zmin = openmc.ZPlane(surface_id=10244, z0=-33.85)
surf133_zmax = openmc.ZPlane(surface_id=10245, z0=31.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf135_xmin = openmc.XPlane(surface_id=10246, x0=-0.03)
surf135_xmax = openmc.XPlane(surface_id=10247, x0=0.03)
surf135_ymin = openmc.YPlane(surface_id=10248, y0=-0.025)
surf135_ymax = openmc.YPlane(surface_id=10249, y0=0.025)
surf135_zmin = openmc.ZPlane(surface_id=10250, z0=-33.85)
surf135_zmax = openmc.ZPlane(surface_id=10251, z0=31.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf137_xmin = openmc.XPlane(surface_id=10252, x0=-0.03)
surf137_xmax = openmc.XPlane(surface_id=10253, x0=0.03)
surf137_ymin = openmc.YPlane(surface_id=10254, y0=-0.025)
surf137_ymax = openmc.YPlane(surface_id=10255, y0=0.025)
surf137_zmin = openmc.ZPlane(surface_id=10256, z0=-33.85)
surf137_zmax = openmc.ZPlane(surface_id=10257, z0=31.65)

# Slot for Cd wire: outer fuel plate
# Rectangular parallelepiped (6 planes)
surf139_xmin = openmc.XPlane(surface_id=10258, x0=-0.03)
surf139_xmax = openmc.XPlane(surface_id=10259, x0=0.03)
surf139_ymin = openmc.YPlane(surface_id=10260, y0=-0.025)
surf139_ymax = openmc.YPlane(surface_id=10261, y0=0.025)
surf139_zmin = openmc.ZPlane(surface_id=10262, z0=-38.35)
surf139_zmax = openmc.ZPlane(surface_id=10263, z0=31.65)

surf201 = openmc.ZCylinder(surface_id=201, r=1.1125)

surf206 = openmc.ZCylinder(surface_id=206, r=1.1125)

surf211 = openmc.ZCylinder(surface_id=211, r=1.1125)

surf216 = openmc.ZCylinder(surface_id=216, r=1.1125)

surf221 = openmc.ZCylinder(surface_id=221, r=1.1125)

surf226 = openmc.ZCylinder(surface_id=226, r=1.1125)

surf231 = openmc.ZCylinder(surface_id=231, r=1.1125)

surf236 = openmc.ZCylinder(surface_id=236, r=1.1125)

surf241 = openmc.ZCylinder(surface_id=241, r=1.1125)

surf246 = openmc.ZCylinder(surface_id=246, r=1.1125)

surf251 = openmc.ZCylinder(surface_id=251, r=1.1125)

surf256 = openmc.ZCylinder(surface_id=256, r=1.1125)

# 1st (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf301_xmin = openmc.XPlane(surface_id=10264, x0=-3.55)
surf301_xmax = openmc.XPlane(surface_id=10265, x0=3.55)
surf301_ymin = openmc.YPlane(surface_id=10266, y0=-0.08)
surf301_ymax = openmc.YPlane(surface_id=10267, y0=0.08)
surf301_zmin = openmc.ZPlane(surface_id=10268, z0=-33.85)
surf301_zmax = openmc.ZPlane(surface_id=10269, z0=39.65)

# 2nd (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf302_xmin = openmc.XPlane(surface_id=10270, x0=-3.55)
surf302_xmax = openmc.XPlane(surface_id=10271, x0=3.55)
surf302_ymin = openmc.YPlane(surface_id=10272, y0=-0.08)
surf302_ymax = openmc.YPlane(surface_id=10273, y0=0.08)
surf302_zmin = openmc.ZPlane(surface_id=10274, z0=-33.85)
surf302_zmax = openmc.ZPlane(surface_id=10275, z0=39.65)

# 3rd (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf303_xmin = openmc.XPlane(surface_id=10276, x0=-3.55)
surf303_xmax = openmc.XPlane(surface_id=10277, x0=3.55)
surf303_ymin = openmc.YPlane(surface_id=10278, y0=-0.08)
surf303_ymax = openmc.YPlane(surface_id=10279, y0=0.08)
surf303_zmin = openmc.ZPlane(surface_id=10280, z0=-33.85)
surf303_zmax = openmc.ZPlane(surface_id=10281, z0=39.65)

# 4th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf304_xmin = openmc.XPlane(surface_id=10282, x0=-3.55)
surf304_xmax = openmc.XPlane(surface_id=10283, x0=3.55)
surf304_ymin = openmc.YPlane(surface_id=10284, y0=-0.08)
surf304_ymax = openmc.YPlane(surface_id=10285, y0=0.08)
surf304_zmin = openmc.ZPlane(surface_id=10286, z0=-33.85)
surf304_zmax = openmc.ZPlane(surface_id=10287, z0=39.65)

# 5th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf305_xmin = openmc.XPlane(surface_id=10288, x0=-3.55)
surf305_xmax = openmc.XPlane(surface_id=10289, x0=3.55)
surf305_ymin = openmc.YPlane(surface_id=10290, y0=-0.08)
surf305_ymax = openmc.YPlane(surface_id=10291, y0=0.08)
surf305_zmin = openmc.ZPlane(surface_id=10292, z0=-33.85)
surf305_zmax = openmc.ZPlane(surface_id=10293, z0=39.65)

# 6th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf306_xmin = openmc.XPlane(surface_id=10294, x0=-3.55)
surf306_xmax = openmc.XPlane(surface_id=10295, x0=3.55)
surf306_ymin = openmc.YPlane(surface_id=10296, y0=-0.08)
surf306_ymax = openmc.YPlane(surface_id=10297, y0=0.08)
surf306_zmin = openmc.ZPlane(surface_id=10298, z0=-33.85)
surf306_zmax = openmc.ZPlane(surface_id=10299, z0=39.65)

# 7th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf307_xmin = openmc.XPlane(surface_id=10300, x0=-3.55)
surf307_xmax = openmc.XPlane(surface_id=10301, x0=3.55)
surf307_ymin = openmc.YPlane(surface_id=10302, y0=-0.08)
surf307_ymax = openmc.YPlane(surface_id=10303, y0=0.08)
surf307_zmin = openmc.ZPlane(surface_id=10304, z0=-33.85)
surf307_zmax = openmc.ZPlane(surface_id=10305, z0=39.65)

# 8th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf308_xmin = openmc.XPlane(surface_id=10306, x0=-3.55)
surf308_xmax = openmc.XPlane(surface_id=10307, x0=3.55)
surf308_ymin = openmc.YPlane(surface_id=10308, y0=-0.08)
surf308_ymax = openmc.YPlane(surface_id=10309, y0=0.08)
surf308_zmin = openmc.ZPlane(surface_id=10310, z0=-33.85)
surf308_zmax = openmc.ZPlane(surface_id=10311, z0=39.65)

# 9th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf309_xmin = openmc.XPlane(surface_id=10312, x0=-3.55)
surf309_xmax = openmc.XPlane(surface_id=10313, x0=3.55)
surf309_ymin = openmc.YPlane(surface_id=10314, y0=-0.08)
surf309_ymax = openmc.YPlane(surface_id=10315, y0=0.08)
surf309_zmin = openmc.ZPlane(surface_id=10316, z0=-33.85)
surf309_zmax = openmc.ZPlane(surface_id=10317, z0=39.65)

# 10th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf310_xmin = openmc.XPlane(surface_id=10318, x0=-3.55)
surf310_xmax = openmc.XPlane(surface_id=10319, x0=3.55)
surf310_ymin = openmc.YPlane(surface_id=10320, y0=-0.08)
surf310_ymax = openmc.YPlane(surface_id=10321, y0=0.08)
surf310_zmin = openmc.ZPlane(surface_id=10322, z0=-33.85)
surf310_zmax = openmc.ZPlane(surface_id=10323, z0=39.65)

# 11th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf311_xmin = openmc.XPlane(surface_id=10324, x0=-3.55)
surf311_xmax = openmc.XPlane(surface_id=10325, x0=3.55)
surf311_ymin = openmc.YPlane(surface_id=10326, y0=-0.08)
surf311_ymax = openmc.YPlane(surface_id=10327, y0=0.08)
surf311_zmin = openmc.ZPlane(surface_id=10328, z0=-33.85)
surf311_zmax = openmc.ZPlane(surface_id=10329, z0=39.65)

# 12th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf312_xmin = openmc.XPlane(surface_id=10330, x0=-3.55)
surf312_xmax = openmc.XPlane(surface_id=10331, x0=3.55)
surf312_ymin = openmc.YPlane(surface_id=10332, y0=-0.08)
surf312_ymax = openmc.YPlane(surface_id=10333, y0=0.08)
surf312_zmin = openmc.ZPlane(surface_id=10334, z0=-33.85)
surf312_zmax = openmc.ZPlane(surface_id=10335, z0=39.65)

# 13th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf313_xmin = openmc.XPlane(surface_id=10336, x0=-3.55)
surf313_xmax = openmc.XPlane(surface_id=10337, x0=3.55)
surf313_ymin = openmc.YPlane(surface_id=10338, y0=-0.08)
surf313_ymax = openmc.YPlane(surface_id=10339, y0=0.08)
surf313_zmin = openmc.ZPlane(surface_id=10340, z0=-33.85)
surf313_zmax = openmc.ZPlane(surface_id=10341, z0=39.65)

# 14th (inner) plate and slot
# Rectangular parallelepiped (6 planes)
surf314_xmin = openmc.XPlane(surface_id=10342, x0=-3.55)
surf314_xmax = openmc.XPlane(surface_id=10343, x0=3.55)
surf314_ymin = openmc.YPlane(surface_id=10344, y0=-0.08)
surf314_ymax = openmc.YPlane(surface_id=10345, y0=0.08)
surf314_zmin = openmc.ZPlane(surface_id=10346, z0=-33.85)
surf314_zmax = openmc.ZPlane(surface_id=10347, z0=39.65)

# Slot for Cd wire: outer fuel plate
# Rectangular parallelepiped (6 planes)
surf321_xmin = openmc.XPlane(surface_id=10348, x0=-0.03)
surf321_xmax = openmc.XPlane(surface_id=10349, x0=0.03)
surf321_ymin = openmc.YPlane(surface_id=10350, y0=-0.025)
surf321_ymax = openmc.YPlane(surface_id=10351, y0=0.025)
surf321_zmin = openmc.ZPlane(surface_id=10352, z0=-38.35)
surf321_zmax = openmc.ZPlane(surface_id=10353, z0=39.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf323_xmin = openmc.XPlane(surface_id=10354, x0=-0.03)
surf323_xmax = openmc.XPlane(surface_id=10355, x0=0.03)
surf323_ymin = openmc.YPlane(surface_id=10356, y0=-0.025)
surf323_ymax = openmc.YPlane(surface_id=10357, y0=0.025)
surf323_zmin = openmc.ZPlane(surface_id=10358, z0=-33.85)
surf323_zmax = openmc.ZPlane(surface_id=10359, z0=39.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf325_xmin = openmc.XPlane(surface_id=10360, x0=-0.03)
surf325_xmax = openmc.XPlane(surface_id=10361, x0=0.03)
surf325_ymin = openmc.YPlane(surface_id=10362, y0=-0.025)
surf325_ymax = openmc.YPlane(surface_id=10363, y0=0.025)
surf325_zmin = openmc.ZPlane(surface_id=10364, z0=-33.85)
surf325_zmax = openmc.ZPlane(surface_id=10365, z0=39.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf327_xmin = openmc.XPlane(surface_id=10366, x0=-0.03)
surf327_xmax = openmc.XPlane(surface_id=10367, x0=0.03)
surf327_ymin = openmc.YPlane(surface_id=10368, y0=-0.025)
surf327_ymax = openmc.YPlane(surface_id=10369, y0=0.025)
surf327_zmin = openmc.ZPlane(surface_id=10370, z0=-33.85)
surf327_zmax = openmc.ZPlane(surface_id=10371, z0=39.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf329_xmin = openmc.XPlane(surface_id=10372, x0=-0.03)
surf329_xmax = openmc.XPlane(surface_id=10373, x0=0.03)
surf329_ymin = openmc.YPlane(surface_id=10374, y0=-0.025)
surf329_ymax = openmc.YPlane(surface_id=10375, y0=0.025)
surf329_zmin = openmc.ZPlane(surface_id=10376, z0=-33.85)
surf329_zmax = openmc.ZPlane(surface_id=10377, z0=39.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf331_xmin = openmc.XPlane(surface_id=10378, x0=-0.03)
surf331_xmax = openmc.XPlane(surface_id=10379, x0=0.03)
surf331_ymin = openmc.YPlane(surface_id=10380, y0=-0.025)
surf331_ymax = openmc.YPlane(surface_id=10381, y0=0.025)
surf331_zmin = openmc.ZPlane(surface_id=10382, z0=-33.85)
surf331_zmax = openmc.ZPlane(surface_id=10383, z0=39.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf333_xmin = openmc.XPlane(surface_id=10384, x0=-0.03)
surf333_xmax = openmc.XPlane(surface_id=10385, x0=0.03)
surf333_ymin = openmc.YPlane(surface_id=10386, y0=-0.025)
surf333_ymax = openmc.YPlane(surface_id=10387, y0=0.025)
surf333_zmin = openmc.ZPlane(surface_id=10388, z0=-33.85)
surf333_zmax = openmc.ZPlane(surface_id=10389, z0=39.65)

# Slot for Cd wire: inner fuel plate
# Rectangular parallelepiped (6 planes)
surf335_xmin = openmc.XPlane(surface_id=10390, x0=-0.03)
surf335_xmax = openmc.XPlane(surface_id=10391, x0=0.03)
surf335_ymin = openmc.YPlane(surface_id=10392, y0=-0.025)
surf335_ymax = openmc.YPlane(surface_id=10393, y0=0.025)
surf335_zmin = openmc.ZPlane(surface_id=10394, z0=-33.85)
surf335_zmax = openmc.ZPlane(surface_id=10395, z0=39.65)

# Al6061 external guide plate slot
# Rectangular parallelepiped (6 planes)
surf341_xmin = openmc.XPlane(surface_id=10396, x0=-3.555)
surf341_xmax = openmc.XPlane(surface_id=10397, x0=3.555)
surf341_ymin = openmc.YPlane(surface_id=10398, y0=-0.065)
surf341_ymax = openmc.YPlane(surface_id=10399, y0=0.065)
surf341_zmin = openmc.ZPlane(surface_id=10400, z0=-38.35)
surf341_zmax = openmc.ZPlane(surface_id=10401, z0=39.65)

# Al6061 internal guide plate slot
# Rectangular parallelepiped (6 planes)
surf342_xmin = openmc.XPlane(surface_id=10402, x0=-3.555)
surf342_xmax = openmc.XPlane(surface_id=10403, x0=3.555)
surf342_ymin = openmc.YPlane(surface_id=10404, y0=-0.065)
surf342_ymax = openmc.YPlane(surface_id=10405, y0=0.065)
surf342_zmin = openmc.ZPlane(surface_id=10406, z0=-33.85)
surf342_zmax = openmc.ZPlane(surface_id=10407, z0=39.65)

# Al6061 internal guide plate slot
# Rectangular parallelepiped (6 planes)
surf343_xmin = openmc.XPlane(surface_id=10408, x0=-3.555)
surf343_xmax = openmc.XPlane(surface_id=10409, x0=3.555)
surf343_ymin = openmc.YPlane(surface_id=10410, y0=-0.065)
surf343_ymax = openmc.YPlane(surface_id=10411, y0=0.065)
surf343_zmin = openmc.ZPlane(surface_id=10412, z0=-33.85)
surf343_zmax = openmc.ZPlane(surface_id=10413, z0=39.65)

# Al6061 external guide plate slot
# Rectangular parallelepiped (6 planes)
surf344_xmin = openmc.XPlane(surface_id=10414, x0=-3.555)
surf344_xmax = openmc.XPlane(surface_id=10415, x0=3.555)
surf344_ymin = openmc.YPlane(surface_id=10416, y0=-0.065)
surf344_ymax = openmc.YPlane(surface_id=10417, y0=0.065)
surf344_zmin = openmc.ZPlane(surface_id=10418, z0=-38.35)
surf344_zmax = openmc.ZPlane(surface_id=10419, z0=39.65)

# Control element, fully extracted: upper
# Rectangular parallelepiped (6 planes)
surf401_xmin = openmc.XPlane(surface_id=10420, x0=-3.245)
surf401_xmax = openmc.XPlane(surface_id=10421, x0=3.245)
surf401_ymin = openmc.YPlane(surface_id=10422, y0=-0.215)
surf401_ymax = openmc.YPlane(surface_id=10423, y0=0.215)
surf401_zmin = openmc.ZPlane(surface_id=10424, z0=30.45)
surf401_zmax = openmc.ZPlane(surface_id=10425, z0=117.15)

# Control element, fully extracted: lower
# Rectangular parallelepiped (6 planes)
surf402_xmin = openmc.XPlane(surface_id=10426, x0=-3.245)
surf402_xmax = openmc.XPlane(surface_id=10427, x0=3.245)
surf402_ymin = openmc.YPlane(surface_id=10428, y0=-0.215)
surf402_ymax = openmc.YPlane(surface_id=10429, y0=0.215)
surf402_zmin = openmc.ZPlane(surface_id=10430, z0=30.45)
surf402_zmax = openmc.ZPlane(surface_id=10431, z0=117.15)

# Control element 4: upper
# Rectangular parallelepiped (6 planes)
surf403_xmin = openmc.XPlane(surface_id=10432, x0=-3.245)
surf403_xmax = openmc.XPlane(surface_id=10433, x0=3.245)
surf403_ymin = openmc.YPlane(surface_id=10434, y0=-0.215)
surf403_ymax = openmc.YPlane(surface_id=10435, y0=0.215)
surf403_zmin = openmc.ZPlane(surface_id=10436, z0=30.45)
surf403_zmax = openmc.ZPlane(surface_id=10437, z0=117.15)

# Control element 4: lower
# Rectangular parallelepiped (6 planes)
surf404_xmin = openmc.XPlane(surface_id=10438, x0=-3.245)
surf404_xmax = openmc.XPlane(surface_id=10439, x0=3.245)
surf404_ymin = openmc.YPlane(surface_id=10440, y0=-0.215)
surf404_ymax = openmc.YPlane(surface_id=10441, y0=0.215)
surf404_zmin = openmc.ZPlane(surface_id=10442, z0=30.45)
surf404_zmax = openmc.ZPlane(surface_id=10443, z0=117.15)

# Rectangular parallelepiped (6 planes)
surf901_xmin = openmc.XPlane(surface_id=10444, x0=-15.4)
surf901_xmax = openmc.XPlane(surface_id=10445, x0=23.1)
surf901_ymin = openmc.YPlane(surface_id=10446, y0=-24.3)
surf901_ymax = openmc.YPlane(surface_id=10447, y0=16.2)
surf901_zmin = openmc.ZPlane(surface_id=10448, z0=-180)
surf901_zmax = openmc.ZPlane(surface_id=10449, z0=800)

surf999 = openmc.ZCylinder(surface_id=999, r=120)


# ==============================================================================
# Universes (from COG define unit blocks)
# ==============================================================================

# Unit 1: Inner fuel plate and slot
u1_cell0 = openmc.Cell(fill=mat1, name="Al-U3Si2")
u1_cell0.region = (+surf21_xmin & -surf21_xmax & +surf21_ymin & -surf21_ymax & +surf21_zmin & -surf21_zmax)
u1_cell1 = openmc.Cell(fill=mat2, name="Al6061")
u1_cell1.region = (-surf21_xmin | +surf21_xmax | -surf21_ymin | +surf21_ymax | -surf21_zmin | +surf21_zmax) & (+surf22_xmin & -surf22_xmax & +surf22_ymin & -surf22_ymax & +surf22_zmin & -surf22_zmax)
u1_cell2 = openmc.Cell(fill=mat3, name="H2O")
u1_cell2.region = (-surf22_xmin | +surf22_xmax | -surf22_ymin | +surf22_ymax | -surf22_zmin | +surf22_zmax) & (+surf23_xmin & -surf23_xmax & +surf23_ymin & -surf23_ymax & +surf23_zmin & -surf23_zmax)
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2])

# Unit 2: Outer fuel plate and slot
u2_cell0 = openmc.Cell(fill=mat1, name="Al-U3Si2")
u2_cell0.region = (+surf21_xmin & -surf21_xmax & +surf21_ymin & -surf21_ymax & +surf21_zmin & -surf21_zmax)
u2_cell1 = openmc.Cell(fill=mat2, name="Al6061")
u2_cell1.region = (-surf21_xmin | +surf21_xmax | -surf21_ymin | +surf21_ymax | -surf21_zmin | +surf21_zmax) & (+surf24_xmin & -surf24_xmax & +surf24_ymin & -surf24_ymax & +surf24_zmin & -surf24_zmax)
u2_cell2 = openmc.Cell(fill=mat3, name="H2O")
u2_cell2.region = (-surf24_xmin | +surf24_xmax | -surf24_ymin | +surf24_ymax | -surf24_zmin | +surf24_zmax) & (+surf25_xmin & -surf25_xmax & +surf25_ymin & -surf25_ymax & +surf25_zmin & -surf25_zmax)
universe2 = openmc.Universe(universe_id=2, cells=[u2_cell0, u2_cell1, u2_cell2])

# Unit 3: Cd wire and slot
u3_cell0 = openmc.Cell(fill=mat4, name="Cd")
u3_cell0.region = -surf26
universe3 = openmc.Universe(universe_id=3, cells=[u3_cell0])

# Unit 4: Normal fuel element
u4_cell0 = openmc.Cell(fill=mat2, name="Al6061")
u4_cell0.region = (-surf41_xmin | +surf41_xmax | -surf41_ymin | +surf41_ymax | -surf41_zmin | +surf41_zmax) & (+surf42_xmin & -surf42_xmax & +surf42_ymin & -surf42_ymax & +surf42_zmin & -surf42_zmax) & (-surf23_xmin | +surf23_xmax | -surf23_ymin | +surf23_ymax | -surf23_zmin | +surf23_zmax) & (-surf101_xmin | +surf101_xmax | -surf101_ymin | +surf101_ymax | -surf101_zmin | +surf101_zmax) & (-surf102_xmin | +surf102_xmax | -surf102_ymin | +surf102_ymax | -surf102_zmin | +surf102_zmax) & (-surf103_xmin | +surf103_xmax | -surf103_ymin | +surf103_ymax | -surf103_zmin | +surf103_zmax) & (-surf104_xmin | +surf104_xmax | -surf104_ymin | +surf104_ymax | -surf104_zmin | +surf104_zmax) & (-surf105_xmin | +surf105_xmax | -surf105_ymin | +surf105_ymax | -surf105_zmin | +surf105_zmax) & (-surf106_xmin | +surf106_xmax | -surf106_ymin | +surf106_ymax | -surf106_zmin | +surf106_zmax) & (-surf107_xmin | +surf107_xmax | -surf107_ymin | +surf107_ymax | -surf107_zmin | +surf107_zmax) & (-surf108_xmin | +surf108_xmax | -surf108_ymin | +surf108_ymax | -surf108_zmin | +surf108_zmax) & (-surf109_xmin | +surf109_xmax | -surf109_ymin | +surf109_ymax | -surf109_zmin | +surf109_zmax)
u4_cell1 = openmc.Cell(fill=mat2, name="Al6061")
u4_cell1.region = (-surf41_xmin | +surf41_xmax | -surf41_ymin | +surf41_ymax | -surf41_zmin | +surf41_zmax) & (-surf42_xmin | +surf42_xmax | -surf42_ymin | +surf42_ymax | -surf42_zmin | +surf42_zmax) & +surf30 & -surf31
universe4 = openmc.Universe(universe_id=4, cells=[u4_cell0, u4_cell1])

# Unit 5: Control fuel element
u5_cell0 = openmc.Cell(fill=mat2, name="Al6061")
u5_cell0.region = (-surf41_xmin | +surf41_xmax | -surf41_ymin | +surf41_ymax | -surf41_zmin | +surf41_zmax) & (+surf42_xmin & -surf42_xmax & +surf42_ymin & -surf42_ymax & +surf42_zmin & -surf42_zmax) & (-surf301_xmin | +surf301_xmax | -surf301_ymin | +surf301_ymax | -surf301_zmin | +surf301_zmax) & (-surf302_xmin | +surf302_xmax | -surf302_ymin | +surf302_ymax | -surf302_zmin | +surf302_zmax) & (-surf303_xmin | +surf303_xmax | -surf303_ymin | +surf303_ymax | -surf303_zmin | +surf303_zmax) & (-surf304_xmin | +surf304_xmax | -surf304_ymin | +surf304_ymax | -surf304_zmin | +surf304_zmax) & (-surf305_xmin | +surf305_xmax | -surf305_ymin | +surf305_ymax | -surf305_zmin | +surf305_zmax) & (-surf306_xmin | +surf306_xmax | -surf306_ymin | +surf306_ymax | -surf306_zmin | +surf306_zmax) & (-surf307_xmin | +surf307_xmax | -surf307_ymin | +surf307_ymax | -surf307_zmin | +surf307_zmax) & (-surf308_xmin | +surf308_xmax | -surf308_ymin | +surf308_ymax | -surf308_zmin | +surf308_zmax) & (-surf309_xmin | +surf309_xmax | -surf309_ymin | +surf309_ymax | -surf309_zmin | +surf309_zmax) & (-surf310_xmin | +surf310_xmax | -surf310_ymin | +surf310_ymax | -surf310_zmin | +surf310_zmax)
u5_cell1 = openmc.Cell(fill=mat2, name="Al6061")
u5_cell1.region = (-surf41_xmin | +surf41_xmax | -surf41_ymin | +surf41_ymax | -surf41_zmin | +surf41_zmax) & (-surf42_xmin | +surf42_xmax | -surf42_ymin | +surf42_ymax | -surf42_zmin | +surf42_zmax) & +surf30 & -surf31
universe5 = openmc.Universe(universe_id=5, cells=[u5_cell0, u5_cell1])

# Unit 6: Grid plate unit cell with no nozzle
u6_cell0 = openmc.Cell(fill=mat5, name="Al")
u6_cell0.region = (+surf90_xmin & -surf90_xmax & +surf90_ymin & -surf90_ymax & +surf90_zmin & -surf90_zmax) & +surf91
u6_cell1 = openmc.Cell(fill=mat3, name="H2O")
u6_cell1.region = (+surf90_xmin & -surf90_xmax & +surf90_ymin & -surf90_ymax & +surf90_zmin & -surf90_zmax) & -surf91
universe6 = openmc.Universe(universe_id=6, cells=[u6_cell0, u6_cell1])

# Unit 7: Grid plate unit cell with a nozzle
u7_cell0 = openmc.Cell(fill=mat5, name="Al")
u7_cell0.region = (+surf90_xmin & -surf90_xmax & +surf90_ymin & -surf90_ymax & +surf90_zmin & -surf90_zmax) & +surf91
universe7 = openmc.Universe(universe_id=7, cells=[u7_cell0])

# Unit 8: Grid plate
# Lattice for unit 8: 8x10 array
# Pitch: (7.700000, 8.100000) cm
# Lower left: (-30.8, -40.5)

# Creating 8x10 lattice with mixed universes
lattice8 = openmc.RectLattice(lattice_id=8)
lattice8.lower_left = [-30.8, -40.5]
lattice8.pitch = [7.700000, 8.100000]
# TODO: Set up universe array for mixed lattice
lattice8.universes = [[universe6]*8]*10
universe8 = lattice8  # Lattice can be used as universe

# Unit 9: Grid plate with primary and secondary holes
universe9 = openmc.Universe(universe_id=9, cells=[])

# Unit 10: Al6061 internal guide plate in water
u10_cell0 = openmc.Cell(fill=mat2, name="Al6061")
u10_cell0.region = (+surf43_xmin & -surf43_xmax & +surf43_ymin & -surf43_ymax & +surf43_zmin & -surf43_zmax)
universe10 = openmc.Universe(universe_id=10, cells=[u10_cell0])

# Unit 11: Al6061 External guide plate in water
u11_cell0 = openmc.Cell(fill=mat2, name="Al6061")
u11_cell0.region = (+surf44_xmin & -surf44_xmax & +surf44_ymin & -surf44_ymax & +surf44_zmin & -surf44_zmax)
universe11 = openmc.Universe(universe_id=11, cells=[u11_cell0])

# Unit 12: Control element, fully withdrawn
u12_cell0 = openmc.Cell(fill=mat6, name="AgInCd")
u12_cell0.region = (+surf50_xmin & -surf50_xmax & +surf50_ymin & -surf50_ymax & +surf50_zmin & -surf50_zmax) & +surf51 & -surf52 & +surf55 & (+surf56_xmin & -surf56_xmax & +surf56_ymin & -surf56_ymax & +surf56_zmin & -surf56_zmax)
u12_cell1 = openmc.Cell(fill=mat6, name="AgInCd")
u12_cell1.region = (+surf50_xmin & -surf50_xmax & +surf50_ymin & -surf50_ymax & +surf50_zmin & -surf50_zmax) & -surf51 & -surf53 & (+surf56_xmin & -surf56_xmax & +surf56_ymin & -surf56_ymax & +surf56_zmin & -surf56_zmax)
u12_cell2 = openmc.Cell(fill=mat6, name="AgInCd")
u12_cell2.region = (+surf50_xmin & -surf50_xmax & +surf50_ymin & -surf50_ymax & +surf50_zmin & -surf50_zmax) & +surf52 & -surf54 & (+surf56_xmin & -surf56_xmax & +surf56_ymin & -surf56_ymax & +surf56_zmin & -surf56_zmax)
u12_cell3 = openmc.Cell(fill=mat7, name="SS304L")
u12_cell3.region = (-surf56_xmin | +surf56_xmax | -surf56_ymin | +surf56_ymax | -surf56_zmin | +surf56_zmax) & (+surf57_xmin & -surf57_xmax & +surf57_ymin & -surf57_ymax & +surf57_zmin & -surf57_zmax)
u12_cell4 = openmc.Cell(fill=mat3, name="H2O")
u12_cell4.region = (-surf57_xmin | +surf57_xmax | -surf57_ymin | +surf57_ymax | -surf57_zmin | +surf57_zmax) & -surf999
universe12 = openmc.Universe(universe_id=12, cells=[u12_cell0, u12_cell1, u12_cell2, u12_cell3, u12_cell4])

# Unit 13: Control fuel elements 1, 2 and 3: fully extracted
universe13 = openmc.Universe(universe_id=13, cells=[])

# Unit 14: Control fuel element 4: partially extracted
universe14 = openmc.Universe(universe_id=14, cells=[])

# Unit 15: Water
u15_cell0 = openmc.Cell(fill=mat3, name="H2O")
u15_cell0.region = -surf999
universe15 = openmc.Universe(universe_id=15, cells=[u15_cell0])

# Unit 16: Core: array of NFEs and CFEs
# Lattice for unit 16: 5x5 array
# Pitch: (7.700000, 8.100000) cm
# Lower left: (-15.4, -24.3)

# Creating 5x5 lattice with mixed universes
lattice16 = openmc.RectLattice(lattice_id=16)
lattice16.lower_left = [-15.4, -24.3]
lattice16.pitch = [7.700000, 8.100000]
# TODO: Set up universe array for mixed lattice
lattice16.universes = [[universe15]*5]*5
universe16 = lattice16  # Lattice can be used as universe

# Cell: H2O
cell0 = openmc.Cell(cell_id=0, fill=mat3, name="H2O")
cell0.region = -surf999 & (-surf901_xmin | +surf901_xmax | -surf901_ymin | +surf901_ymax | -surf901_zmin | +surf901_zmax) & (-surf90_xmin | +surf90_xmax | -surf90_ymin | +surf90_ymax | -surf90_zmin | +surf90_zmax) & (-surf60_xmin | +surf60_xmax | -surf60_ymin | +surf60_ymax | -surf60_zmin | +surf60_zmax)

# Cell: H2O
cell1 = openmc.Cell(cell_id=1, fill=mat3, name="H2O")
cell1.region = (-surf22_xmin | +surf22_xmax | -surf22_ymin | +surf22_ymax | -surf22_zmin | +surf22_zmax) & (+surf23_xmin & -surf23_xmax & +surf23_ymin & -surf23_ymax & +surf23_zmin & -surf23_zmax)

# Cell: H2O
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="H2O")
cell2.region = (-surf24_xmin | +surf24_xmax | -surf24_ymin | +surf24_ymax | -surf24_zmin | +surf24_zmax) & (+surf25_xmin & -surf25_xmax & +surf25_ymin & -surf25_ymax & +surf25_zmin & -surf25_zmax)

# Cell: H2O
cell3 = openmc.Cell(cell_id=3, fill=mat3, name="H2O")
cell3.region = (+surf90_xmin & -surf90_xmax & +surf90_ymin & -surf90_ymax & +surf90_zmin & -surf90_zmax) & -surf91

# Cell: H2O
cell4 = openmc.Cell(cell_id=4, fill=mat3, name="H2O")
cell4.region = (-surf57_xmin | +surf57_xmax | -surf57_ymin | +surf57_ymax | -surf57_zmin | +surf57_zmax) & -surf999

# Cell: H2O
cell5 = openmc.Cell(cell_id=5, fill=mat3, name="H2O")
cell5.region = -surf999

# Cell: Cd
cell6 = openmc.Cell(cell_id=6, fill=mat4, name="Cd")
cell6.region = -surf999 & (+surf60_xmin & -surf60_xmax & +surf60_ymin & -surf60_ymax & +surf60_zmin & -surf60_zmax) & -surf61 & +surf62

# Cell: Cd
cell7 = openmc.Cell(cell_id=7, fill=mat4, name="Cd")
cell7.region = -surf999 & (+surf60_xmin & -surf60_xmax & +surf60_ymin & -surf60_ymax & +surf60_zmin & -surf60_zmax) & -surf63 & +surf64

# Cell: Al
cell8 = openmc.Cell(cell_id=8, fill=mat5, name="Al")
cell8.region = -surf999 & (+surf60_xmin & -surf60_xmax & +surf60_ymin & -surf60_ymax & +surf60_zmin & -surf60_zmax) & +surf61

# Cell: Al
cell9 = openmc.Cell(cell_id=9, fill=mat5, name="Al")
cell9.region = -surf999 & (+surf60_xmin & -surf60_xmax & +surf60_ymin & -surf60_ymax & +surf60_zmin & -surf60_zmax) & -surf62 & +surf63

# Cell: Al
cell10 = openmc.Cell(cell_id=10, fill=mat5, name="Al")
cell10.region = (+surf90_xmin & -surf90_xmax & +surf90_ymin & -surf90_ymax & +surf90_zmin & -surf90_zmax) & +surf91

# Cell: Al2O3
cell11 = openmc.Cell(cell_id=11, fill=mat8, name="Al2O3")
cell11.region = -surf999 & (+surf60_xmin & -surf60_xmax & +surf60_ymin & -surf60_ymax & +surf60_zmin & -surf60_zmax) & -surf64

# Cell using unit 9: GP
cell12 = openmc.Cell(cell_id=12, fill=universe9, name="GP")
cell12.region = -surf999 & (+surf90_xmin & -surf90_xmax & +surf90_ymin & -surf90_ymax & +surf90_zmin & -surf90_zmax) & (-surf60_xmin | +surf60_xmax | -surf60_ymin | +surf60_ymax | -surf60_zmin | +surf60_zmax)

# Cell using unit 16: Core
cell13 = openmc.Cell(cell_id=13, fill=universe16, name="Core")
cell13.region = -surf999 & (+surf901_xmin & -surf901_xmax & +surf901_ymin & -surf901_ymax & +surf901_zmin & -surf901_zmax) & (-surf90_xmin | +surf90_xmax | -surf90_ymin | +surf90_ymax | -surf90_zmin | +surf90_zmax)

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary (6 planes)
# TODO: Adjust dimensions to encompass your entire geometry
boundary_xmin = openmc.XPlane(surface_id=10450, x0=-200, boundary_type="vacuum")
boundary_xmax = openmc.XPlane(surface_id=10451, x0=200, boundary_type="vacuum")
boundary_ymin = openmc.YPlane(surface_id=10452, y0=-200, boundary_type="vacuum")
boundary_ymax = openmc.YPlane(surface_id=10453, y0=200, boundary_type="vacuum")
boundary_zmin = openmc.ZPlane(surface_id=10454, z0=-200, boundary_type="vacuum")
boundary_zmax = openmc.ZPlane(surface_id=10455, z0=200, boundary_type="vacuum")

# Create outer void cell (everything outside geometry but inside boundary)
# Particles are killed at the vacuum boundary
outer_region = +boundary_xmin & -boundary_xmax & +boundary_ymin & -boundary_ymax & +boundary_zmin & -boundary_zmax
outer_region = outer_region & ~cell0.region
outer_region = outer_region & ~cell1.region
outer_region = outer_region & ~cell2.region
outer_region = outer_region & ~cell3.region
outer_region = outer_region & ~cell4.region
outer_region = outer_region & ~cell5.region
outer_region = outer_region & ~cell6.region
outer_region = outer_region & ~cell7.region
outer_region = outer_region & ~cell8.region
outer_region = outer_region & ~cell9.region
outer_region = outer_region & ~cell10.region
outer_region = outer_region & ~cell11.region
outer_region = outer_region & ~cell12.region
outer_region = outer_region & ~cell13.region
outer_cell = openmc.Cell(cell_id=14, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, outer_cell])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 1000
settings.batches = 5050
settings.inactive = 51
settings.run_mode = "eigenvalue"

# Source definition
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)
settings.source = source

# Enable delayed neutron kinetics and alpha eigenvalue calculations
settings.calculate_prompt_k = True
settings.calculate_alpha = True

settings.export_to_xml()

# ==============================================================================
# Tallies
# ==============================================================================

tallies = openmc.Tallies()
tallies.export_to_xml()

# ==============================================================================
# Run OpenMC
# ==============================================================================

openmc.run()
