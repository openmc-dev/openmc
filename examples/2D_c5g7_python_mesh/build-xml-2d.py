import openmc
import openmc.mgxs
import numpy as np

z_min = -32.13
z_max = 32.13

###############################################################################
# General Settings
###############################################################################

# Instantiate a Settings, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.energy_mode = "multi-group"
settings_file.batches = 20 #1000
settings_file.inactive = 10 #600
settings_file.particles = 1750
settings_file.run_mode = 'eigenvalue'
settings_file.output = {'tallies': False, 'summary': False}

# Create an initial uniform spatial source distribution for sampling rays.
# Note that this must be uniform in space and angle.
lower_left = (-64.26/2, -64.26/2, -1e-5)
upper_right = (64.26/2, 64.26/2, 1e-5)
uniform_dist = openmc.stats.Box(lower_left, upper_right)
settings_file.random_ray['ray_source'] = openmc.IndependentSource(space=uniform_dist)
settings_file.random_ray['distance_inactive'] = 60.0
settings_file.random_ray['distance_active'] = 250.0
settings_file.random_ray['source_shape'] = 'flat'
settings_file.random_ray['sample_method'] = 'prng'

# A mesh with pincell level resolution will be overlaid by default.
# Increasing this parameter above 1 will subdivide the mesh further.
mesh_resolution_multiplier = 7


###############################################################################
# Plots
###############################################################################

plot_1 = openmc.Plot(plot_id=1)
plot_1.filename = 'plot_1'
plot_1.origin = [0.0, 0.0, 0.0]
plot_1.width = [64.26, 64.26, 1.0]
plot_1.pixels = [1000, 1000, 1]
# plot_1.pixels = [3000, 3000, 1]
plot_1.type = 'voxel'

# Instantiate a Plots collection and export to XML
plot_file = openmc.Plots([plot_1])

###############################################################################
# Pinmaker Utility
###############################################################################

def pinmaker(inner_fill, outer_fill, fuel_radius):

    pincell_base = openmc.Universe()

    cylinder = openmc.ZCylinder(r=fuel_radius)
    inside_cell = openmc.Cell(fill=inner_fill, region=-cylinder)
    outside_cell = openmc.Cell(fill=outer_fill, region=+cylinder)
    pincell_base.add_cells([inside_cell, outside_cell])

    return pincell_base

###############################################################################
# MGXS
###############################################################################

# Instantiate the energy group data
groups = openmc.mgxs.EnergyGroups(np.logspace(-5, 7, 8))

# Instantiate the 7-group C5G7 cross section data
uo2_xsdata = openmc.XSdata('uo2', groups)
uo2_xsdata.order = 0
uo2_xsdata.set_total(
    np.array([1.779490E-01, 3.298050E-01, 4.803880E-01, 5.543670E-01,
              3.118010E-01, 3.951680E-01, 5.644060E-01]))
uo2_xsdata.set_absorption(
    np.array([8.02480E-03, 3.71740E-03, 2.67690E-02, 9.62360E-02, 3.00200E-02,
              1.11260E-01, 2.82780E-01]))
scatter_matrix = \
    [[[1.275370E-01, 4.237800E-02, 9.437400E-06, 5.516300E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 3.244560E-01, 1.631400E-03, 3.142700E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 4.509400E-01, 2.679200E-03, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 4.525650E-01, 5.566400E-03, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 1.252500E-04, 2.714010E-01, 1.025500E-02, 1.002100E-08],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 1.296800E-03, 2.658020E-01, 1.680900E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 8.545800E-03, 2.730800E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
uo2_xsdata.set_scatter_matrix(scatter_matrix)
uo2_xsdata.set_fission(
    np.array([7.212060E-03, 8.193010E-04, 6.453200E-03, 1.856480E-02, 1.780840E-02, 8.303480E-02, 2.160040E-01]))
uo2_xsdata.set_nu_fission(
    np.array([2.005998E-02, 2.027303E-03, 1.570599E-02, 4.518301E-02, 4.334208E-02, 2.020901E-01, 5.257105E-01]))
uo2_xsdata.set_chi(
    np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]))

mox43_xsdata = openmc.XSdata('mox43', groups)
mox43_xsdata.order = 0
mox43_xsdata.set_total(
    np.array([1.787310E-01, 3.308490E-01, 4.837720E-01, 5.669220E-01, 4.262270E-01, 6.789970E-01, 6.828520E-01]))
mox43_xsdata.set_absorption(
    np.array([8.43390E-03, 3.75770E-03, 2.79700E-02, 1.04210E-01, 1.39940E-01, 4.09180E-01, 4.09350E-01]))
scatter_matrix = \
    [[[1.288760E-01, 4.141300E-02, 8.229000E-06, 5.040500E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 3.254520E-01, 1.639500E-03, 1.598200E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 4.531880E-01, 2.614200E-03, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 4.571730E-01, 5.539400E-03, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 1.604600E-04, 2.768140E-01, 9.312700E-03, 9.165600E-09],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 2.005100E-03, 2.529620E-01, 1.485000E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 8.494800E-03, 2.650070E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
mox43_xsdata.set_scatter_matrix(scatter_matrix)
mox43_xsdata.set_fission(
    np.array([7.62704E-03, 8.76898E-04, 5.69835E-03, 2.28872E-02, 1.07635E-02, 2.32757E-01, 2.48968E-01]))
mox43_xsdata.set_nu_fission(
    np.array([2.175300E-02, 2.535103E-03, 1.626799E-02, 6.547410E-02, 3.072409E-02, 6.666510E-01, 7.139904E-01]))
mox43_xsdata.set_chi(
    np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]))

mox7_xsdata = openmc.XSdata('mox7', groups)
mox7_xsdata.order = 0
mox7_xsdata.set_total(
    np.array([1.813230E-01, 3.343680E-01, 4.937850E-01, 5.912160E-01, 4.741980E-01, 8.336010E-01, 8.536030E-01]))
mox7_xsdata.set_absorption(
    np.array([9.06570E-03, 4.29670E-03, 3.28810E-02, 1.22030E-01, 1.82980E-01, 5.68460E-01, 5.85210E-01]))
scatter_matrix = \
    [[[1.304570E-01, 4.179200E-02, 8.510500E-06, 5.132900E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 3.284280E-01, 1.643600E-03, 2.201700E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 4.583710E-01, 2.533100E-03, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 4.637090E-01, 5.476600E-03, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 1.761900E-04, 2.823130E-01, 8.728900E-03, 9.001600E-09],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 2.276000E-03, 2.497510E-01, 1.311400E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 8.864500E-03, 2.595290E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
mox7_xsdata.set_scatter_matrix(scatter_matrix)
mox7_xsdata.set_fission(
    np.array([8.25446E-03, 1.32565E-03, 8.42156E-03, 3.28730E-02, 1.59636E-02, 3.23794E-01, 3.62803E-01]))
mox7_xsdata.set_nu_fission(
    np.array([2.381395E-02, 3.858689E-03, 2.413400E-02, 9.436622E-02, 4.576988E-02, 9.281814E-01, 1.043200E+00]))
mox7_xsdata.set_chi(
    np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]))

mox87_xsdata = openmc.XSdata('mox87', groups)
mox87_xsdata.order = 0
mox87_xsdata.set_total(
    np.array([1.830450E-01, 3.367050E-01, 5.005070E-01, 6.061740E-01, 5.027540E-01, 9.210280E-01, 9.552310E-01]))
mox87_xsdata.set_absorption(
    np.array([9.48620E-03, 4.65560E-03, 3.62400E-02, 1.32720E-01, 2.08400E-01, 6.58700E-01, 6.90170E-01]))
scatter_matrix = \
    [[[1.315040E-01, 4.204600E-02, 8.697200E-06, 5.193800E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 3.304030E-01, 1.646300E-03, 2.600600E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 4.617920E-01, 2.474900E-03, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 4.680210E-01, 5.433000E-03, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 1.859700E-04, 2.857710E-01, 8.397300E-03, 8.928000E-09],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 2.391600E-03, 2.476140E-01, 1.232200E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 8.968100E-03, 2.560930E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
mox87_xsdata.set_scatter_matrix(scatter_matrix)
mox87_xsdata.set_fission(
    np.array([8.67209E-03, 1.62426E-03, 1.02716E-02, 3.90447E-02, 1.92576E-02, 3.74888E-01, 4.30599E-01]))
mox87_xsdata.set_nu_fission(
    np.array([2.518600E-02, 4.739509E-03, 2.947805E-02, 1.122500E-01, 5.530301E-02, 1.074999E+00, 1.239298E+00]))
mox87_xsdata.set_chi(
    np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]))

fiss_chamber_xsdata = openmc.XSdata('fiss_chamber', groups)
fiss_chamber_xsdata.order = 0
fiss_chamber_xsdata.set_total(
    np.array([1.260320E-01, 2.931600E-01, 2.842500E-01, 2.810200E-01, 3.344600E-01, 5.656400E-01, 1.172140E+00]))
fiss_chamber_xsdata.set_absorption(
    np.array([5.11320E-04, 7.58130E-05, 3.16430E-04, 1.16750E-03, 3.39770E-03, 9.18860E-03, 2.32440E-02]))
scatter_matrix = \
    [[[6.616590E-02, 5.907000E-02, 2.833400E-04, 1.462200E-06, 2.064200E-08, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 2.403770E-01, 5.243500E-02, 2.499000E-04, 1.923900E-05, 2.987500E-06, 4.214000E-07],
      [0.000000E-00, 0.000000E-00, 1.834250E-01, 9.228800E-02, 6.936500E-03, 1.079000E-03, 2.054300E-04],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 7.907690E-02, 1.699900E-01, 2.586000E-02, 4.925600E-03],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 3.734000E-05, 9.975700E-02, 2.067900E-01, 2.447800E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 9.174200E-04, 3.167740E-01, 2.387600E-01],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 4.979300E-02, 1.09910E+00]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
fiss_chamber_xsdata.set_scatter_matrix(scatter_matrix)
fiss_chamber_xsdata.set_fission(
    np.array([4.79002E-09, 5.82564E-09, 4.63719E-07, 5.24406E-06, 1.45390E-07, 7.14972E-07, 2.08041E-06]))
fiss_chamber_xsdata.set_nu_fission(
    np.array([1.323401E-08, 1.434500E-08, 1.128599E-06, 1.276299E-05, 3.538502E-07, 1.740099E-06, 5.063302E-06]))
fiss_chamber_xsdata.set_chi(
    np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]))

guide_tube_xsdata = openmc.XSdata('guide_tube', groups)
guide_tube_xsdata.order = 0
guide_tube_xsdata.set_total(
    np.array([1.260320E-01, 2.931600E-01, 2.842400E-01, 2.809600E-01, 3.344400E-01, 5.656400E-01, 1.172150E+00]))
guide_tube_xsdata.set_absorption(
    np.array([5.11320E-04, 7.58010E-05, 3.15720E-04, 1.15820E-03, 3.39750E-03, 9.18780E-03, 2.32420E-02]))
scatter_matrix = \
    [[[6.616590E-02, 5.907000E-02, 2.833400E-04, 1.462200E-06, 2.064200E-08, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 2.403770E-01, 5.243500E-02, 2.499000E-04, 1.923900E-05, 2.987500E-06, 4.214000E-07],
      [0.000000E-00, 0.000000E-00, 1.832970E-01, 9.239700E-02, 6.944600E-03, 1.080300E-03, 2.056700E-04],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 7.885110E-02, 1.701400E-01, 2.588100E-02, 4.929700E-03],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 3.733300E-05, 9.973720E-02, 2.067900E-01, 2.447800E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 9.172600E-04, 3.167650E-01, 2.387700E-01],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 4.979200E-02, 1.099120E+00]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
guide_tube_xsdata.set_scatter_matrix(scatter_matrix)

water_xsdata = openmc.XSdata('water', groups)
water_xsdata.order = 0
water_xsdata.set_total(
    np.array([1.592060E-01, 4.129700E-01, 5.903100E-01, 5.843500E-01, 7.180000E-01, 1.254450E+00, 2.650380E+00]))
water_xsdata.set_absorption(
    np.array([6.01050E-04, 1.57930E-05, 3.37160E-04, 1.94060E-03, 5.74160E-03, 1.50010E-02, 3.72390E-02]))
scatter_matrix = \
    [[[4.447770E-02, 1.134000E-01, 7.234700E-04, 3.749900E-06, 5.318400E-08, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 2.823340E-01, 1.299400E-01, 6.234000E-04, 4.800200E-05, 7.448600E-06, 1.045500E-06],
      [0.000000E-00, 0.000000E-00, 3.452560E-01, 2.245700E-01, 1.699900E-02, 2.644300E-03, 5.034400E-04],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 9.102840E-02, 4.155100E-01, 6.373200E-02, 1.213900E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 7.143700E-05, 1.391380E-01, 5.118200E-01, 6.122900E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 2.215700E-03, 6.999130E-01, 5.373200E-01],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 1.324400E-01, 2.480700E+00]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
water_xsdata.set_scatter_matrix(scatter_matrix)

control_rod_xsdata = openmc.XSdata('control_rod', groups)
control_rod_xsdata.order = 0
control_rod_xsdata.set_total(
    np.array([2.16768E-01, 4.80098E-01, 8.86369E-01, 9.70009E-01, 9.10482E-01, 1.13775E+00, 1.84048E+00]))
control_rod_xsdata.set_absorption(
    np.array([1.70490E-03, 8.36224E-03, 8.37901E-02, 3.97797E-01, 6.98763E-01, 9.29508E-01, 1.17836E+00]))
scatter_matrix = \
    [[[1.70563E-01, 4.44012E-02, 9.83670E-05, 1.27786E-07, 0.00000E-00, 0.00000E-00, 0.00000E-00],
      [0.00000E-00, 4.71050E-01, 6.85480E-04, 3.91395E-10, 0.00000E-00, 0.00000E-00, 0.00000E-00],
      [0.00000E-00, 0.00000E-00, 8.01859E-01, 7.20132E-04, 0.00000E-00, 0.00000E-00, 0.00000E-00],
      [0.00000E-00, 0.00000E-00, 0.00000E-00, 5.70752E-01, 1.46015E-03, 0.00000E-00, 0.00000E-00],
      [0.00000E-00, 0.00000E-00, 0.00000E-00, 6.55562E-05, 2.07838E-01, 3.81486E-03, 3.69760E-09],
      [0.00000E-00, 0.00000E-00, 0.00000E-00, 0.00000E-00, 1.02427E-03, 2.02465E-01, 4.75290E-03],
      [0.00000E-00, 0.00000E-00, 0.00000E-00, 0.00000E-00, 0.00000E-00, 3.53043E-03, 6.58597E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
control_rod_xsdata.set_scatter_matrix(scatter_matrix)

mg_cross_sections_file = openmc.MGXSLibrary(groups)
mg_cross_sections_file.add_xsdatas([uo2_xsdata, mox43_xsdata, mox7_xsdata, mox87_xsdata,
                                    fiss_chamber_xsdata, guide_tube_xsdata, water_xsdata,
                                    control_rod_xsdata])
mg_cross_sections_file.export_to_hdf5()

###############################################################################
# Materials
###############################################################################

# Instantiate some Macroscopic Data
uo2_data = openmc.Macroscopic('uo2')
mox43_data = openmc.Macroscopic('mox43')
mox7_data = openmc.Macroscopic('mox7')
mox87_data = openmc.Macroscopic('mox87')
fiss_chamber_data = openmc.Macroscopic('fiss_chamber')
guide_tube_data = openmc.Macroscopic('guide_tube')
water_data = openmc.Macroscopic('water')
control_rod_data = openmc.Macroscopic('control_rod')

# Instantiate Materials dictionary
materials = {}

# Instantiate some Materials and register the appropriate Nuclides
materials['UO2'] = openmc.Material(name='UO2')
materials['UO2'].set_density('macro', 1.0)
materials['UO2'].add_macroscopic(uo2_data)

materials['MOX 4.3%'] = openmc.Material(name='MOX 4.3%')
materials['MOX 4.3%'].set_density('macro', 1.0)
materials['MOX 4.3%'].add_macroscopic(mox43_data)

materials['MOX 7.0%'] = openmc.Material(name='MOX 7.0%')
materials['MOX 7.0%'].set_density('macro', 1.0)
materials['MOX 7.0%'].add_macroscopic(mox7_data)

materials['MOX 8.7%'] = openmc.Material(name='MOX 8.7%')
materials['MOX 8.7%'].set_density('macro', 1.0)
materials['MOX 8.7%'].add_macroscopic(mox87_data)

materials['Fission Chamber'] = openmc.Material(name='Fission Chamber')
materials['Fission Chamber'].set_density('macro', 1.0)
materials['Fission Chamber'].add_macroscopic(fiss_chamber_data)

materials['Guide Tube'] = openmc.Material(name='Guide Tube')
materials['Guide Tube'].set_density('macro', 1.0)
materials['Guide Tube'].add_macroscopic(guide_tube_data)

materials['Water'] = openmc.Material(name='Water')
materials['Water'].set_density('macro', 1.0)
materials['Water'].add_macroscopic(water_data)

materials['Control Rod'] = openmc.Material(name='Control Rod')
materials['Control Rod'].set_density('macro', 1.0)
materials['Control Rod'].add_macroscopic(control_rod_data)

# Instantiate a Materials collection, register all Materials, and export to XML
materials_file = openmc.Materials(materials.values())
materials_file.cross_sections = './mgxs.h5'

###############################################################################
# Surfaces
###############################################################################

# Create a dictionary to store the surfaces
surfaces = {}

# Instantiate Pin Cell ZCylinder surface
surfaces['Pin Cell ZCylinder'] = openmc.ZCylinder(x0=0, y0=0, r=0.54, name='Pin Cell ZCylinder')
surfaces['Core x-min']         = openmc.XPlane(x0=-32.13, name='Core x-min')
surfaces['Core x-max']         = openmc.XPlane(x0= 32.13, name='Core x-max')
surfaces['Core y-min']         = openmc.YPlane(y0=-32.13, name='Core y-min')
surfaces['Core y-max']         = openmc.YPlane(y0= 32.13, name='Core y-max')
surfaces['Small Core z-min']   = openmc.ZPlane(z0=z_min, name='Small Core z-min')
surfaces['Small Core z-max']   = openmc.ZPlane(z0=z_max, name='Small Core z-max')
surfaces['Big Core z-min']     = openmc.ZPlane(z0=-107.1, name='Big Core z-min')
surfaces['Big Core z-max']     = openmc.ZPlane(z0= 107.1, name='Big Core z-max')
surfaces['Pin x-min']          = openmc.XPlane(x0=-0.63, name='Pin x-min')
surfaces['Pin x-max']          = openmc.XPlane(x0=+0.63, name='Pin x-max')
surfaces['Pin y-min']          = openmc.YPlane(y0=-0.63, name='Pin y-min')
surfaces['Pin y-max']          = openmc.YPlane(y0=+0.63, name='Pin y-max')

surfaces['Core x-max'].boundary_type       = 'vacuum'
surfaces['Core y-min'].boundary_type       = 'vacuum'
surfaces['Core y-max'].boundary_type       = 'reflective'
surfaces['Small Core z-min'].boundary_type = 'reflective'
surfaces['Small Core z-max'].boundary_type = 'reflective'
surfaces['Big Core z-min'].boundary_type   = 'reflective'
surfaces['Big Core z-max'].boundary_type   = 'vacuum'
surfaces['Pin x-min'].boundary_type        = 'reflective'
surfaces['Pin x-max'].boundary_type        = 'reflective'
surfaces['Pin y-min'].boundary_type        = 'reflective'
surfaces['Pin y-max'].boundary_type        = 'reflective'
surfaces['Core x-min'].boundary_type       = 'reflective'

###############################################################################
# Cells
###############################################################################

# Instantiate Cells
cells = {}
cells['UO2']                         = openmc.Cell(name='UO2')
cells['MOX 4.3%']                    = openmc.Cell(name='MOX 4.3%')
cells['MOX 7.0%']                    = openmc.Cell(name='MOX 7.0%')
cells['MOX 8.7%']                    = openmc.Cell(name='MOX 8.7%')
cells['Fission Chamber']             = openmc.Cell(name='Fission Chamber')
cells['Guide Tube']                  = openmc.Cell(name='Guide Tube')
cells['Reflector']                   = openmc.Cell(name='Reflector')
cells['Control Rod']                 = openmc.Cell(name='Control Rod')
cells['UO2 Moderator']               = openmc.Cell(name='UO2 Moderator')
cells['MOX 4.3% Moderator']          = openmc.Cell(name='MOX 4.3% Moderator')
cells['MOX 7.0% Moderator']          = openmc.Cell(name='MOX 7.0% Moderator')
cells['MOX 8.7% Moderator']          = openmc.Cell(name='MOX 8.7% Moderator')
cells['Fission Chamber Moderator']   = openmc.Cell(name='Fission Chamber Moderator')
cells['Guide Tube Moderator']        = openmc.Cell(name='Guide Tube Moderator')
cells['Control Rod Moderator']       = openmc.Cell(name='Control Rod Moderator')
cells['UO2 Unrodded Assembly']       = openmc.Cell(name='UO2 Unrodded Assembly')
cells['UO2 Rodded Assembly']         = openmc.Cell(name='UO2 Rodded Assembly')
cells['MOX Unrodded Assembly']       = openmc.Cell(name='MOX Unrodded Assembly')
cells['MOX Rodded Assembly']         = openmc.Cell(name='MOX Rodded Assembly')
cells['Reflector Unrodded Assembly'] = openmc.Cell(name='Water Unrodded Assembly')
cells['Reflector Rodded Assembly']   = openmc.Cell(name='Water Rodded Assembly')
cells['Core']                        = openmc.Cell(name='Core')
cells['UO2 Pin']                     = openmc.Cell(name='UO2 Pin')

cells['Reflector'] = openmc.Cell(name='Reflector')

# Use surface half-spaces to define regions
cells['UO2'].region                       = -surfaces['Pin Cell ZCylinder']
cells['MOX 4.3%'].region                  = -surfaces['Pin Cell ZCylinder']
cells['MOX 7.0%'].region                  = -surfaces['Pin Cell ZCylinder']
cells['MOX 8.7%'].region                  = -surfaces['Pin Cell ZCylinder']
cells['Fission Chamber'].region           = -surfaces['Pin Cell ZCylinder']
cells['Guide Tube'].region                = -surfaces['Pin Cell ZCylinder']
cells['Control Rod'].region               = -surfaces['Pin Cell ZCylinder']
cells['UO2 Moderator'].region             = +surfaces['Pin Cell ZCylinder']
cells['MOX 4.3% Moderator'].region        = +surfaces['Pin Cell ZCylinder']
cells['MOX 7.0% Moderator'].region        = +surfaces['Pin Cell ZCylinder']
cells['MOX 8.7% Moderator'].region        = +surfaces['Pin Cell ZCylinder']
cells['Fission Chamber Moderator'].region = +surfaces['Pin Cell ZCylinder']
cells['Guide Tube Moderator'].region      = +surfaces['Pin Cell ZCylinder']
cells['Control Rod Moderator'].region     = +surfaces['Pin Cell ZCylinder']
cells['UO2 Pin'].region                   = +surfaces['Pin Cell ZCylinder'] & \
                                            +surfaces['Pin x-min'] & \
                                            -surfaces['Pin x-max'] & \
                                            +surfaces['Pin y-min'] & \
                                            -surfaces['Pin y-max']

# Register Materials with Cells
cells['UO2'].fill                       = materials['UO2']
cells['MOX 4.3%'].fill                  = materials['MOX 4.3%']
cells['MOX 7.0%'].fill                  = materials['MOX 7.0%']
cells['MOX 8.7%'].fill                  = materials['MOX 8.7%']
cells['Fission Chamber'].fill           = materials['Fission Chamber']
cells['Guide Tube'].fill                = materials['Guide Tube']
cells['Control Rod'].fill               = materials['Control Rod']
cells['UO2 Moderator'].fill             = materials['Water']
cells['MOX 4.3% Moderator'].fill        = materials['Water']
cells['MOX 7.0% Moderator'].fill        = materials['Water']
cells['MOX 8.7% Moderator'].fill        = materials['Water']
cells['Fission Chamber Moderator'].fill = materials['Water']
cells['Guide Tube Moderator'].fill      = materials['Water']
cells['Control Rod Moderator'].fill     = materials['Water']
cells['Reflector'].fill                 = materials['Water']
cells['UO2 Pin'].fill                   = materials['Water']

###############################################################################
# Universes
###############################################################################

# Instantiate Universes
universes = {}
universes['Root']                        = openmc.Universe(universe_id=0,  name='Root')
universes['UO2']                         = openmc.Universe(universe_id=1,  name='UO2')
universes['MOX 4.3%']                    = openmc.Universe(universe_id=2,  name='MOX 4.3%')
universes['MOX 7.0%']                    = openmc.Universe(universe_id=3,  name='MOX 7.0%')
universes['MOX 8.7%']                    = openmc.Universe(universe_id=4,  name='MOX 8.7%')
universes['Fission Chamber']             = openmc.Universe(universe_id=5,  name='Fission Chamber')
universes['Guide Tube']                  = openmc.Universe(universe_id=6,  name='Guide Tube')
universes['Control Rod']                 = openmc.Universe(universe_id=7,  name='Control Rod')
universes['Reflector']                   = openmc.Universe(universe_id=8,  name='Reflector')
universes['UO2 Unrodded Assembly']       = openmc.Universe(universe_id=9,  name='UO2 Unrodded Assembly')
universes['UO2 Rodded Assembly']         = openmc.Universe(universe_id=10, name='UO2 Rodded Assembly')
universes['MOX Unrodded Assembly']       = openmc.Universe(universe_id=11, name='MOX Unrodded Assembly')
universes['MOX Rodded Assembly']         = openmc.Universe(universe_id=12, name='MOX Rodded Assembly')
universes['Reflector Unrodded Assembly'] = openmc.Universe(universe_id=13, name='Reflector Unrodded Assembly')
universes['Reflector Rodded Assembly']   = openmc.Universe(universe_id=14, name='Reflector Rodded Assembly')

universes['Reflector']            = openmc.Universe(name='Reflector')


# Register Cells with Universes
universes['Root']                       .add_cell(cells['Core'])
universes['UO2']                        .add_cells([cells['UO2'], cells['UO2 Moderator']])
universes['MOX 4.3%']                   .add_cells([cells['MOX 4.3%'], cells['MOX 4.3% Moderator']])
universes['MOX 7.0%']                   .add_cells([cells['MOX 7.0%'], cells['MOX 7.0% Moderator']])
universes['MOX 8.7%']                   .add_cells([cells['MOX 8.7%'], cells['MOX 8.7% Moderator']])
universes['Fission Chamber']            .add_cells([cells['Fission Chamber'], cells['Fission Chamber Moderator']])
universes['Guide Tube']                 .add_cells([cells['Guide Tube'], cells['Guide Tube Moderator']])
universes['Control Rod']                .add_cells([cells['Control Rod'], cells['Control Rod Moderator']])
universes['Reflector']                  .add_cell(cells['Reflector'])
universes['UO2 Unrodded Assembly']      .add_cell(cells['UO2 Unrodded Assembly'])
universes['UO2 Rodded Assembly']        .add_cell(cells['UO2 Rodded Assembly'])
universes['MOX Unrodded Assembly']      .add_cell(cells['MOX Unrodded Assembly'])
universes['MOX Rodded Assembly']        .add_cell(cells['MOX Rodded Assembly'])
universes['Reflector Unrodded Assembly'].add_cell(cells['Reflector Unrodded Assembly'])
universes['Reflector Rodded Assembly']  .add_cell(cells['Reflector Rodded Assembly'])

universes['Reflector']           .add_cell(cells['Reflector'])

###############################################################################
# Subdivided Pincells
###############################################################################


pins_to_make = ['UO2', 'MOX 4.3%', 'MOX 7.0%', 'MOX 8.7%', 'Fission Chamber', 'Guide Tube', 'Control Rod']
for pin in pins_to_make:
    universes[pin] = pinmaker(materials[pin], materials['Water'], 0.54)

###############################################################################
#                     Create a dictionary of the assembly lattices
###############################################################################

# Instantiate the Lattices
lattices = {}
lattices['UO2 Unrodded Assembly'] = \
    openmc.RectLattice(lattice_id=101, name='UO2 Unrodded Assembly')
lattices['UO2 Unrodded Assembly'].dimension = [17, 17]
lattices['UO2 Unrodded Assembly'].lower_left = [-10.71, -10.71]
lattices['UO2 Unrodded Assembly'].pitch = [1.26, 1.26]
u = universes['UO2']
g = universes['Guide Tube']
f = universes['Fission Chamber']
lattices['UO2 Unrodded Assembly'].universes = \
    [[u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, g, u, u, g, u, u, g, u, u, u, u, u],
     [u, u, u, g, u, u, u, u, u, u, u, u, u, g, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, g, u, u, g, u, u, g, u, u, g, u, u, g, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, g, u, u, g, u, u, f, u, u, g, u, u, g, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, g, u, u, g, u, u, g, u, u, g, u, u, g, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, g, u, u, u, u, u, u, u, u, u, g, u, u, u],
     [u, u, u, u, u, g, u, u, g, u, u, g, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u]]

lattices['UO2 Rodded Assembly'] = \
    openmc.RectLattice(lattice_id=102, name='UO2 Rodded Assembly')
lattices['UO2 Rodded Assembly'].dimension = [17, 17]
lattices['UO2 Rodded Assembly'].lower_left = [-10.71, -10.71]
lattices['UO2 Rodded Assembly'].pitch = [1.26, 1.26]
u = universes['UO2']
r = universes['Control Rod']
f = universes['Fission Chamber']
lattices['UO2 Rodded Assembly'].universes = \
    [[u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, r, u, u, r, u, u, r, u, u, u, u, u],
     [u, u, u, r, u, u, u, u, u, u, u, u, u, r, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, r, u, u, r, u, u, r, u, u, r, u, u, r, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, r, u, u, r, u, u, f, u, u, r, u, u, r, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, r, u, u, r, u, u, r, u, u, r, u, u, r, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, r, u, u, u, u, u, u, u, u, u, r, u, u, u],
     [u, u, u, u, u, r, u, u, r, u, u, r, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u]]

lattices['MOX Unrodded Assembly'] = \
    openmc.RectLattice(lattice_id=103, name='MOX Unrodded Assembly')
lattices['MOX Unrodded Assembly'].dimension = [17, 17]
lattices['MOX Unrodded Assembly'].lower_left = [-10.71, -10.71]
lattices['MOX Unrodded Assembly'].pitch = [1.26, 1.26]
m = universes['MOX 4.3%']
n = universes['MOX 7.0%']
o = universes['MOX 8.7%']
g = universes['Guide Tube']
f = universes['Fission Chamber']
lattices['MOX Unrodded Assembly'].universes = \
    [[m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m],
     [m, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, m],
     [m, n, n, n, n, g, n, n, g, n, n, g, n, n, n, n, m],
     [m, n, n, g, n, o, o, o, o, o, o, o, n, g, n, n, m],
     [m, n, n, n, o, o, o, o, o, o, o, o, o, n, n, n, m],
     [m, n, g, o, o, g, o, o, g, o, o, g, o, o, g, n, m],
     [m, n, n, o, o, o, o, o, o, o, o, o, o, o, n, n, m],
     [m, n, n, o, o, o, o, o, o, o, o, o, o, o, n, n, m],
     [m, n, g, o, o, g, o, o, f, o, o, g, o, o, g, n, m],
     [m, n, n, o, o, o, o, o, o, o, o, o, o, o, n, n, m],
     [m, n, n, o, o, o, o, o, o, o, o, o, o, o, n, n, m],
     [m, n, g, o, o, g, o, o, g, o, o, g, o, o, g, n, m],
     [m, n, n, n, o, o, o, o, o, o, o, o, o, n, n, n, m],
     [m, n, n, g, n, o, o, o, o, o, o, o, n, g, n, n, m],
     [m, n, n, n, n, g, n, n, g, n, n, g, n, n, n, n, m],
     [m, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, m],
     [m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m]]

lattices['MOX Rodded Assembly'] = \
    openmc.RectLattice(lattice_id=104, name='MOX Rodded Assembly')
lattices['MOX Rodded Assembly'].dimension = [17, 17]
lattices['MOX Rodded Assembly'].lower_left = [-10.71, -10.71]
lattices['MOX Rodded Assembly'].pitch = [1.26, 1.26]
m = universes['MOX 4.3%']
n = universes['MOX 7.0%']
o = universes['MOX 8.7%']
r = universes['Control Rod']
f = universes['Fission Chamber']
lattices['MOX Rodded Assembly'].universes = \
    [[m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m],
     [m, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, m],
     [m, n, n, n, n, r, n, n, r, n, n, r, n, n, n, n, m],
     [m, n, n, r, n, o, o, o, o, o, o, o, n, r, n, n, m],
     [m, n, n, n, o, o, o, o, o, o, o, o, o, n, n, n, m],
     [m, n, r, o, o, r, o, o, r, o, o, r, o, o, r, n, m],
     [m, n, n, o, o, o, o, o, o, o, o, o, o, o, n, n, m],
     [m, n, n, o, o, o, o, o, o, o, o, o, o, o, n, n, m],
     [m, n, r, o, o, r, o, o, f, o, o, r, o, o, r, n, m],
     [m, n, n, o, o, o, o, o, o, o, o, o, o, o, n, n, m],
     [m, n, n, o, o, o, o, o, o, o, o, o, o, o, n, n, m],
     [m, n, r, o, o, r, o, o, r, o, o, r, o, o, r, n, m],
     [m, n, n, n, o, o, o, o, o, o, o, o, o, n, n, n, m],
     [m, n, n, r, n, o, o, o, o, o, o, o, n, r, n, n, m],
     [m, n, n, n, n, r, n, n, r, n, n, r, n, n, n, n, m],
     [m, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, m],
     [m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m]]

lattices['Reflector Unrodded Assembly'] = \
    openmc.RectLattice(lattice_id=105, name='Reflector Unrodded Assembly')
lattices['Reflector Unrodded Assembly'].dimension = [1, 1]
lattices['Reflector Unrodded Assembly'].lower_left = [-10.71, -10.71]
lattices['Reflector Unrodded Assembly'].pitch = [21.42, 21.42]
w = universes['Reflector']
lattices['Reflector Unrodded Assembly'].universes = [[w]]



lattices['Reflector Rodded Assembly'] = \
    openmc.RectLattice(lattice_id=106, name='Reflector Rodded Assembly')
lattices['Reflector Rodded Assembly'].dimension = [17, 17]
lattices['Reflector Rodded Assembly'].lower_left = [-10.71, -10.71]
lattices['Reflector Rodded Assembly'].pitch = [1.26, 1.26]
u = universes['Reflector']
r = universes['Control Rod']
lattices['Reflector Rodded Assembly'].universes = \
    [[u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, r, u, u, r, u, u, r, u, u, u, u, u],
     [u, u, u, r, u, u, u, u, u, u, u, u, u, r, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, r, u, u, r, u, u, r, u, u, r, u, u, r, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, r, u, u, r, u, u, u, u, u, r, u, u, r, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, r, u, u, r, u, u, r, u, u, r, u, u, r, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, r, u, u, u, u, u, u, u, u, u, r, u, u, u],
     [u, u, u, u, u, r, u, u, r, u, u, r, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
     [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u]]

# Add lattice to cells
cells['UO2 Unrodded Assembly'].fill       = lattices['UO2 Unrodded Assembly']
cells['UO2 Rodded Assembly'].fill         = lattices['UO2 Rodded Assembly']
cells['MOX Unrodded Assembly'].fill       = lattices['MOX Unrodded Assembly']
cells['MOX Rodded Assembly'].fill         = lattices['MOX Rodded Assembly']
cells['Reflector Unrodded Assembly'].fill = lattices['Reflector Unrodded Assembly']
cells['Reflector Rodded Assembly'].fill   = lattices['Reflector Rodded Assembly']

###############################################################################
# Base Universe
###############################################################################

# Instantiate Core boundaries
#cells['Core'].region = +surfaces['Core x-min'] & +surfaces['Core y-min'] & \
#    -surfaces['Core x-max'] & -surfaces['Core y-max']

cells['Core'].region = +surfaces['Core x-min'] & +surfaces['Core y-min'] & \
    -surfaces['Core x-max'] & -surfaces['Core y-max'] & +surfaces['Small Core z-min'] & -surfaces['Small Core z-max']

lattices['Core'] = openmc.RectLattice(lattice_id=201, name='3x3 core lattice')
lattices['Core'].dimension = [3, 3]
lattices['Core'].lower_left = [-32.13, -32.13]
lattices['Core'].pitch = [21.42, 21.42]
w = universes['Reflector Unrodded Assembly']
u = universes['UO2 Unrodded Assembly']
m = universes['MOX Unrodded Assembly']
r = universes['Reflector']
lattices['Core'].universes = [[u, m, r], [m, u, r], [r, r, r]]
cells['Core'].fill = lattices['Core']

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry()
geometry.root_universe = universes['Root']
geometry.remove_redundant_surfaces()



###############################################################################
# Tallies
###############################################################################

tallies = {}

# Instantiate a tally mesh
mesh = openmc.RegularMesh(mesh_id=1)
mesh.dimension = [51, 51]
mesh.lower_left = [-32.13, -32.13]
mesh.upper_right = [32.13, 32.13]

# Instantiate some tally Filters
mesh_filter = openmc.MeshFilter(mesh)

# Instantiate the Tally
tallies['Mesh Rates'] = openmc.Tally(name='power')
tallies['Mesh Rates'].filters = [mesh_filter]
tallies['Mesh Rates'].scores = ['fission',]

# Instantiate a Tallies, register Tally/Mesh, and export to XML
tallies_file = openmc.Tallies(tallies.values())

###############################################################################
# Random Ray Mesh Overlay for SR Subdivision
###############################################################################

subdivision_mesh = openmc.RegularMesh()
subdivision_mesh.dimension = [d * mesh_resolution_multiplier for d in mesh.dimension]
subdivision_mesh.lower_left = mesh.lower_left
subdivision_mesh.upper_right = mesh.upper_right

domain = geometry.root_universe
settings_file.random_ray['source_region_meshes'] = [(subdivision_mesh, [domain])]

###############################################################################
# Model
###############################################################################

model = openmc.model.Model()
model.geometry = geometry
model.materials = materials_file
model.settings = settings_file
model.xs_data = mg_cross_sections_file
model.tallies = tallies_file
model.plots = plot_file

model.export_to_model_xml()
