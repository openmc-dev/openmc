"""An example file showing how to make a geometry.

This particular example creates a 3x3 geometry, with 8 regular pins and one
Gd-157 2 wt-percent enriched.  All pins are segmented.
"""

from collections import OrderedDict
import math

import numpy as np
import openmc


def density_to_mat(dens_dict):
    """Generates an OpenMC material from a cell ID and self.number_density.

    Parameters
    ----------
    dens_dict : dict
        Dictionary mapping nuclide names to densities

    Returns
    -------
    openmc.Material
        The OpenMC material filled with nuclides.

    """
    mat = openmc.Material()
    for key in dens_dict:
        mat.add_nuclide(key, 1.0e-24*dens_dict[key])
    mat.set_density('sum')

    return mat


def generate_initial_number_density():
    """ Generates initial number density.

    These results were from a CASMO5 run in which the gadolinium pin was
    loaded with 2 wt percent of Gd-157.
    """

    # Concentration to be used for all fuel pins
    fuel_dict = OrderedDict()
    fuel_dict['U235'] = 1.05692e21
    fuel_dict['U234'] = 1.00506e19
    fuel_dict['U238'] = 2.21371e22
    fuel_dict['O16'] = 4.62954e22
    fuel_dict['O17'] = 1.127684e20
    fuel_dict['Xe135'] = 1.0e10
    fuel_dict['Xe136'] = 1.0e10
    fuel_dict['Gd156'] = 1.0e10
    fuel_dict['Gd157'] = 1.0e10
    # fuel_dict['O18'] = 9.51352e19 # Does not exist in ENDF71, merged into 17

    # Concentration to be used for the gadolinium fuel pin
    fuel_gd_dict = OrderedDict()
    fuel_gd_dict['U235'] = 1.03579e21
    fuel_gd_dict['U238'] = 2.16943e22
    fuel_gd_dict['Gd156'] = 3.95517E+10
    fuel_gd_dict['Gd157'] = 1.08156e20
    fuel_gd_dict['O16'] = 4.64035e22
    fuel_dict['Xe136'] = 1.0e10
    fuel_dict['Xe135'] = 1.0e10
    # There are a whole bunch of 1e-10 stuff here.

    # Concentration to be used for cladding
    clad_dict = OrderedDict()
    clad_dict['O16'] = 3.07427e20
    clad_dict['O17'] = 7.48868e17
    clad_dict['Cr50'] = 3.29620e18
    clad_dict['Cr52'] = 6.35639e19
    clad_dict['Cr53'] = 7.20763e18
    clad_dict['Cr54'] = 1.79413e18
    clad_dict['Fe54'] = 5.57350e18
    clad_dict['Fe56'] = 8.74921e19
    clad_dict['Fe57'] = 2.02057e18
    clad_dict['Fe58'] = 2.68901e17
    clad_dict['Cr50'] = 3.29620e18
    clad_dict['Cr52'] = 6.35639e19
    clad_dict['Cr53'] = 7.20763e18
    clad_dict['Cr54'] = 1.79413e18
    clad_dict['Ni58'] = 2.51631e19
    clad_dict['Ni60'] = 9.69278e18
    clad_dict['Ni61'] = 4.21338e17
    clad_dict['Ni62'] = 1.34341e18
    clad_dict['Ni64'] = 3.43127e17
    clad_dict['Zr90'] = 2.18320e22
    clad_dict['Zr91'] = 4.76104e21
    clad_dict['Zr92'] = 7.27734e21
    clad_dict['Zr94'] = 7.37494e21
    clad_dict['Zr96'] = 1.18814e21
    clad_dict['Sn112'] = 4.67352e18
    clad_dict['Sn114'] = 3.17992e18
    clad_dict['Sn115'] = 1.63814e18
    clad_dict['Sn116'] = 7.00546e19
    clad_dict['Sn117'] = 3.70027e19
    clad_dict['Sn118'] = 1.16694e20
    clad_dict['Sn119'] = 4.13872e19
    clad_dict['Sn120'] = 1.56973e20
    clad_dict['Sn122'] = 2.23076e19
    clad_dict['Sn124'] = 2.78966e19

    # Gap concentration
    # Funny enough, the example problem uses air.
    gap_dict = OrderedDict()
    gap_dict['O16'] = 7.86548e18
    gap_dict['O17'] = 2.99548e15
    gap_dict['N14'] = 3.38646e19
    gap_dict['N15'] = 1.23717e17

    # Concentration to be used for coolant
    # No boron
    cool_dict = OrderedDict()
    cool_dict['H1'] = 4.68063e22
    cool_dict['O16'] = 2.33427e22
    cool_dict['O17'] = 8.89086e18

    # Store these dictionaries in the initial conditions dictionary
    initial_density = OrderedDict()
    initial_density['fuel_gd'] = fuel_gd_dict
    initial_density['fuel'] = fuel_dict
    initial_density['gap'] = gap_dict
    initial_density['clad'] = clad_dict
    initial_density['cool'] = cool_dict

    # Set up libraries to use
    temperature = OrderedDict()
    sab = OrderedDict()

    # Toggle betweeen MCNP and NNDC data
    MCNP = False

    if MCNP:
        temperature['fuel_gd'] = 900.0
        temperature['fuel'] = 900.0
        # We approximate temperature of everything as 600K, even though it was
        # actually 580K.
        temperature['gap'] = 600.0
        temperature['clad'] = 600.0
        temperature['cool'] = 600.0
    else:
        temperature['fuel_gd'] = 293.6
        temperature['fuel'] = 293.6
        temperature['gap'] = 293.6
        temperature['clad'] = 293.6
        temperature['cool'] = 293.6

    sab['cool'] = 'c_H_in_H2O'

    # Set up burnable materials
    burn = OrderedDict()
    burn['fuel_gd'] = True
    burn['fuel'] = True
    burn['gap'] = False
    burn['clad'] = False
    burn['cool'] = False

    return temperature, sab, initial_density, burn


def segment_pin(n_rings, n_wedges, r_fuel, r_gap, r_clad):
    """ Calculates a segmented pin.

    Separates a pin with n_rings and n_wedges.  All cells have equal volume.
    Pin is centered at origin.
    """

    # Calculate all the volumes of interest
    v_fuel = math.pi * r_fuel**2
    v_gap = math.pi * r_gap**2 - v_fuel
    v_clad = math.pi * r_clad**2 - v_fuel - v_gap
    v_ring = v_fuel / n_rings
    v_segment = v_ring / n_wedges

    # Compute ring radiuses
    r_rings = np.zeros(n_rings)

    for i in range(n_rings):
        r_rings[i] = math.sqrt(1.0/(math.pi) * v_ring * (i+1))

    # Compute thetas
    theta = np.linspace(0, 2*math.pi, n_wedges + 1)

    # Compute surfaces
    fuel_rings = [openmc.ZCylinder(x0=0, y0=0, r=r_rings[i])
                  for i in range(n_rings)]

    fuel_wedges = [openmc.Plane(a=math.cos(theta[i]), b=math.sin(theta[i]))
                   for i in range(n_wedges)]

    gap_ring = openmc.ZCylinder(x0=0, y0=0, r=r_gap)
    clad_ring = openmc.ZCylinder(x0=0, y0=0, r=r_clad)

    # Create cells
    fuel_cells = []
    if n_wedges == 1:
        for i in range(n_rings):
            cell = openmc.Cell(name='fuel')
            if i == 0:
                cell.region = -fuel_rings[0]
            else:
                cell.region = +fuel_rings[i-1] & -fuel_rings[i]
            fuel_cells.append(cell)
    else:
        for i in range(n_rings):
            for j in range(n_wedges):
                cell = openmc.Cell(name='fuel')
                if i == 0:
                    if j != n_wedges-1:
                        cell.region = (-fuel_rings[0]
                                       & +fuel_wedges[j]
                                       & -fuel_wedges[j+1])
                    else:
                        cell.region = (-fuel_rings[0]
                                       & +fuel_wedges[j]
                                       & -fuel_wedges[0])
                else:
                    if j != n_wedges-1:
                        cell.region = (+fuel_rings[i-1]
                                       & -fuel_rings[i]
                                       & +fuel_wedges[j]
                                       & -fuel_wedges[j+1])
                    else:
                        cell.region = (+fuel_rings[i-1]
                                       & -fuel_rings[i]
                                       & +fuel_wedges[j]
                                       & -fuel_wedges[0])
                fuel_cells.append(cell)

    # Gap ring
    gap_cell = openmc.Cell(name='gap')
    gap_cell.region = +fuel_rings[-1] & -gap_ring
    fuel_cells.append(gap_cell)

    # Clad ring
    clad_cell = openmc.Cell(name='clad')
    clad_cell.region = +gap_ring & -clad_ring
    fuel_cells.append(clad_cell)

    # Moderator
    mod_cell = openmc.Cell(name='cool')
    mod_cell.region = +clad_ring
    fuel_cells.append(mod_cell)

    # Form universe
    fuel_u = openmc.Universe()
    fuel_u.add_cells(fuel_cells)

    return fuel_u, v_segment, v_gap, v_clad


def generate_geometry(n_rings, n_wedges):
    """ Generates example geometry.

    This function creates the initial geometry, a 9 pin reflective problem.
    One pin, containing gadolinium, is discretized into sectors.

    In addition to what one would do with the general OpenMC geometry code, it
    is necessary to create a dictionary, volume, that maps a cell ID to a
    volume. Further, by naming cells the same as the above materials, the code
    can automatically handle the mapping.

    Parameters
    ----------
    n_rings : int
        Number of rings to generate for the geometry
    n_wedges : int
        Number of wedges to generate for the geometry
    """

    pitch = 1.26197
    r_fuel = 0.412275
    r_gap = 0.418987
    r_clad = 0.476121

    n_pin = 3

    # This table describes the 'fuel' to actual type mapping
    # It's not necessary to do it this way.  Just adjust the initial conditions
    # below.
    mapping = ['fuel', 'fuel', 'fuel',
               'fuel', 'fuel_gd', 'fuel',
               'fuel', 'fuel', 'fuel']

    # Form pin cell
    fuel_u, v_segment, v_gap, v_clad = segment_pin(n_rings, n_wedges, r_fuel, r_gap, r_clad)

    # Form lattice
    all_water_c = openmc.Cell(name='cool')
    all_water_u = openmc.Universe(cells=(all_water_c, ))

    lattice = openmc.RectLattice()
    lattice.pitch = [pitch]*2
    lattice.lower_left = [-pitch*n_pin/2, -pitch*n_pin/2]
    lattice_array = [[fuel_u for i in range(n_pin)] for j in range(n_pin)]
    lattice.universes = lattice_array
    lattice.outer = all_water_u

    # Bound universe
    x_low = openmc.XPlane(-pitch*n_pin/2, boundary_type='reflective')
    x_high = openmc.XPlane(pitch*n_pin/2, boundary_type='reflective')
    y_low = openmc.YPlane(-pitch*n_pin/2, boundary_type='reflective')
    y_high = openmc.YPlane(pitch*n_pin/2, boundary_type='reflective')
    z_low = openmc.ZPlane(-10, boundary_type='reflective')
    z_high = openmc.ZPlane(10, boundary_type='reflective')

    # Compute bounding box
    lower_left = [-pitch*n_pin/2, -pitch*n_pin/2, -10]
    upper_right = [pitch*n_pin/2, pitch*n_pin/2, 10]

    root_c = openmc.Cell(fill=lattice)
    root_c.region = (+x_low & -x_high
                     & +y_low & -y_high
                     & +z_low & -z_high)
    root_u = openmc.Universe(universe_id=0, cells=(root_c, ))
    geometry = openmc.Geometry(root_u)

    v_cool = pitch**2 - (v_gap + v_clad + n_rings * n_wedges * v_segment)

    # Store volumes for later usage
    volume = {'fuel': v_segment, 'gap': v_gap, 'clad': v_clad, 'cool': v_cool}

    return geometry, volume, mapping, lower_left, upper_right


def generate_problem(n_rings=5, n_wedges=8):
    """ Merges geometry and materials.

    This function initializes the materials for each cell using the dictionaries
    provided by generate_initial_number_density.  It is assumed a cell named
    'fuel' will have further region differentiation (see mapping).

    Parameters
    ----------
    n_rings : int, optional
        Number of rings to generate for the geometry
    n_wedges : int, optional
        Number of wedges to generate for the geometry
    """

    # Get materials dictionary, geometry, and volumes
    temperature, sab, initial_density, burn = generate_initial_number_density()
    geometry, volume, mapping, lower_left, upper_right = generate_geometry(n_rings, n_wedges)

    # Apply distribmats, fill geometry
    cells = geometry.root_universe.get_all_cells()
    for cell_id in cells:
        cell = cells[cell_id]
        if cell.name == 'fuel':

            omc_mats = []

            for cell_type in mapping:
                omc_mat = density_to_mat(initial_density[cell_type])

                if cell_type in sab:
                    omc_mat.add_s_alpha_beta(sab[cell_type])
                omc_mat.temperature = temperature[cell_type]
                omc_mat.depletable = burn[cell_type]
                omc_mat.volume = volume['fuel']

                omc_mats.append(omc_mat)

            cell.fill = omc_mats
        elif cell.name != '':
            omc_mat = density_to_mat(initial_density[cell.name])

            if cell.name in sab:
                omc_mat.add_s_alpha_beta(sab[cell.name])
            omc_mat.temperature = temperature[cell.name]
            omc_mat.depletable = burn[cell.name]
            omc_mat.volume = volume[cell.name]

            cell.fill = omc_mat

    return geometry, lower_left, upper_right
