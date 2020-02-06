#!/usr/bin/env python3
"""Update OpenMC's deprecated multi-group cross section XML files to the latest
HDF5-based format.

"""

import os
import warnings
import xml.etree.ElementTree as ET

import argparse
import numpy as np

import openmc.mgxs_library


def parse_args():
    """Read the input files from the commandline."""
    # Create argument parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='input XML file')
    parser.add_argument('-o', '--output', nargs='?', default='',
                        help='output file, in HDF5 format')
    args = vars(parser.parse_args())

    if args['output'] == '':
        filename = args['input'].name
        extension = os.path.splitext(filename)
        if extension == '.xml':
            filename = filename[:filename.rfind('.')] + '.h5'
        args['output'] = filename

    # Parse and return commandline arguments.
    return args


def get_data(element, entry):
    value = element.find(entry)
    if value is not None:
        value = value.text.strip()
    elif entry in element.attrib:
        value = element.attrib[entry].strip()
    else:
        value = None

    return value


def main():
    args = parse_args()

    # Parse the XML data.
    tree = ET.parse(args['input'])
    root = tree.getroot()

    # Get old metadata
    group_structure = tree.find('group_structure').text.strip()
    group_structure = np.array(group_structure.split(), dtype=float)
    # Convert from MeV to eV
    group_structure *= 1.e6
    energy_groups = openmc.mgxs.EnergyGroups(group_structure)

    inverse_velocity = tree.find('inverse-velocity')
    if inverse_velocity is not None:
        inverse_velocity = inverse_velocity.text.split()
        inverse_velocity = np.array(inverse_velocity, dtype=float)
    else:
        inverse_velocity = None

    xsd = []
    names = []

    # Now move on to the cross section data itself
    for xsdata_elem in root.iter('xsdata'):
        name = get_data(xsdata_elem, 'name')

        temperature = get_data(xsdata_elem, 'kT')
        if temperature is not None:
            temperature = float(temperature) / openmc.data.K_BOLTZMANN * 1.E6
        else:
            temperature = 294.
        temperatures = [temperature]

        awr = get_data(xsdata_elem, 'awr')
        if awr is not None:
            awr = float(awr)

        representation = get_data(xsdata_elem, 'representation')
        if representation is None:
            representation = 'isotropic'
        if representation == 'angle':
            n_azi = int(get_data(xsdata_elem, 'num_azimuthal'))
            n_pol = int(get_data(xsdata_elem, 'num_polar'))

        scatter_format = get_data(xsdata_elem, 'scatt_type')
        if scatter_format is None:
            scatter_format = 'legendre'

        order = int(get_data(xsdata_elem, 'order'))

        tab_leg = get_data(xsdata_elem, 'tabular_legendre')
        if tab_leg is not None:
            warnings.warn('The tabular_legendre option has moved to the '
                          'settings.xml file and must be added manually')

        # Either add the data to a previously existing xsdata (if it is
        # for the same 'name' but a different temperature), or create a
        # new one.
        try:
            # It is in our list, so store that entry
            i = names.index(name)
        except ValueError:
            # It is not in our list, so add it
            i = -1
            xsd.append(openmc.XSdata(name, energy_groups,
                                     temperatures=temperatures,
                                     representation=representation))
            if awr is not None:
                xsd[-1].atomic_weight_ratio = awr
            if representation == 'angle':
                xsd[-1].num_azimuthal = n_azi
                xsd[-1].num_polar = n_pol
            xsd[-1].scatter_format = scatter_format
            xsd[-1].order = order
            names.append(name)

        if scatter_format == 'legendre':
            order_dim = order + 1
        else:
            order_dim = order

        if i != -1:
            xsd[i].add_temperature(temperature)

        total = get_data(xsdata_elem, 'total')
        if total is not None:
            total = np.array(total.split(), dtype=float)
            total.shape = xsd[i].xs_shapes['[G]']
            xsd[i].set_total(total, temperature)

        if inverse_velocity is not None:
            xsd[i].set_inverse_velocity(inverse_velocity, temperature)

        absorption = get_data(xsdata_elem, 'absorption')
        absorption = np.array(absorption.split(), dtype=float)
        absorption.shape = xsd[i].xs_shapes['[G]']
        xsd[i].set_absorption(absorption, temperature)

        scatter = get_data(xsdata_elem, 'scatter')
        scatter = np.array(scatter.split(), dtype=float)
        # This is now a flattened-array of something that started with a
        # shape of [Order][G][G']; we need to unflatten and then switch the
        # ordering
        in_shape = (order_dim, energy_groups.num_groups,
                    energy_groups.num_groups)
        if representation == 'angle':
            in_shape = (n_pol, n_azi) + in_shape
            scatter.shape = in_shape
            scatter = np.swapaxes(scatter, 2, 3)
            scatter = np.swapaxes(scatter, 3, 4)
        else:
            scatter.shape = in_shape
            scatter = np.swapaxes(scatter, 0, 1)
            scatter = np.swapaxes(scatter, 1, 2)

        xsd[i].set_scatter_matrix(scatter, temperature)

        multiplicity = get_data(xsdata_elem, 'multiplicity')
        if multiplicity is not None:
            multiplicity = np.array(multiplicity.split(), dtype=float)
            multiplicity.shape = xsd[i].xs_shapes["[G][G']"]
            xsd[i].set_multiplicity_matrix(multiplicity, temperature)

        fission = get_data(xsdata_elem, 'fission')
        if fission is not None:
            fission = np.array(fission.split(), dtype=float)
            fission.shape = xsd[i].xs_shapes['[G]']
            xsd[i].set_fission(fission, temperature)

        kappa_fission = get_data(xsdata_elem, 'kappa_fission')
        if kappa_fission is not None:
            kappa_fission = np.array(kappa_fission.split(), dtype=float)
            kappa_fission.shape = xsd[i].xs_shapes['[G]']
            xsd[i].set_kappa_fission(kappa_fission, temperature)

        chi = get_data(xsdata_elem, 'chi')
        if chi is not None:
            chi = np.array(chi.split(), dtype=float)
            chi.shape = xsd[i].xs_shapes['[G]']
            xsd[i].set_chi(chi, temperature)
        else:
            chi = None

        nu_fission = get_data(xsdata_elem, 'nu_fission')
        if nu_fission is not None:
            nu_fission = np.array(nu_fission.split(), dtype=float)
            if chi is not None:
                nu_fission.shape = xsd[i].xs_shapes['[G]']
            else:
                nu_fission.shape = xsd[i].xs_shapes["[G][G']"]
            xsd[i].set_nu_fission(nu_fission, temperature)

    # Build library as we go, but first we have enough to initialize it
    lib = openmc.MGXSLibrary(energy_groups)
    lib.add_xsdatas(xsd)
    lib.export_to_hdf5(args['output'])


if __name__ == '__main__':
    main()
