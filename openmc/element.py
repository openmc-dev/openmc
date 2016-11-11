from collections import OrderedDict
from numbers import Real
import re
import sys
import os

from six import string_types
from xml.etree import ElementTree as ET
import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.data import NATURAL_ABUNDANCE, atomic_mass
from openmc.plot_data import *


class Element(object):
    """A natural element that auto-expands to add the isotopes of an element to
    a material in their natural abundance. Internally, the OpenMC Python API
    expands the natural element into isotopes only when the materials.xml file
    is created.

    Parameters
    ----------
    name : str
        Chemical symbol of the element, e.g. Pu

    Attributes
    ----------
    name : str
        Chemical symbol of the element, e.g. Pu
    scattering : {'data', 'iso-in-lab', None}
        The type of angular scattering distribution to use

    """

    def __init__(self, name=''):
        # Initialize class attributes
        self._name = ''
        self._scattering = None

        # Set class attributes
        self.name = name

    def __eq__(self, other):
        if isinstance(other, Element):
            if self.name != other.name:
                return False
            else:
                return True
        elif isinstance(other, string_types) and other == self.name:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        return repr(self) > repr(other)

    def __lt__(self, other):
        return not self > other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Element    -    {0}\n'.format(self._name)
        if self.scattering is not None:
            string += '{0: <16}{1}{2}\n'.format('\tscattering', '=\t',
                                                self.scattering)

        return string

    @property
    def name(self):
        return self._name

    @property
    def scattering(self):
        return self._scattering

    @name.setter
    def name(self, name):
        cv.check_type('element name', name, string_types)
        cv.check_length('element name', name, 1, 2)
        self._name = name

    @scattering.setter
    def scattering(self, scattering):

        if not scattering in ['data', 'iso-in-lab', None]:
            msg = 'Unable to set scattering for Element to {0} which ' \
                  'is not "data", "iso-in-lab", or None'.format(scattering)
            raise ValueError(msg)

        self._scattering = scattering

    def expand(self, percent, percent_type, enrichment=None,
               cross_sections=None):
        """Expand natural element into its naturally-occurring isotopes.

        An optional cross_sections argument or the OPENMC_CROSS_SECTIONS
        environment variable is used to specify a cross_sections.xml file.
        If the cross_sections.xml file is found, the element is expanded only
        into the isotopes/nuclides present in cross_sections.xml. If no
        cross_sections.xml file is found, the element is expanded based on its
        naturally occurring isotopes.

        Parameters
        ----------
        percent : float
            Atom or weight percent
        percent_type : {'ao', 'wo'}
            'ao' for atom percent and 'wo' for weight percent
        enrichment : float, optional
            Enrichment for U235 in weight percent. For example, input 4.95 for
            4.95 weight percent enriched U. Default is None
            (natural composition).
        cross_sections : str, optional
            Location of cross_sections.xml file. Default is None.

        Returns
        -------
        isotopes : list
            Naturally-occurring isotopes of the element. Each item of the list
            is a tuple consisting of an openmc.Nuclide instance and the natural
            abundance of the isotope.

        """

        # Get the nuclides present in nature
        natural_nuclides = set()
        for nuclide in sorted(NATURAL_ABUNDANCE.keys()):
            if re.match(r'{}\d+'.format(self.name), nuclide):
                natural_nuclides.add(nuclide)

        # Create dict to store the expanded nuclides and abundances
        abundances = OrderedDict()

        # If cross_sections is None, get the cross sections from the
        # OPENMC_CROSS_SECTIONS environment variable
        if cross_sections is None:
            cross_sections = os.environ.get('OPENMC_CROSS_SECTIONS')

        # If a cross_sections library is present, check natural nuclides
        # against the nuclides in the library
        if cross_sections is not None:

            library_nuclides = set()
            tree = ET.parse(cross_sections)
            root = tree.getroot()
            for child in root:
                nuclide = child.attrib['materials']
                if re.match(r'{}\d+'.format(self.name), nuclide) and \
                   '_m' not in nuclide:
                    library_nuclides.add(nuclide)

            # Get a set of the mutual and absent nuclides. Convert to lists
            # and sort to avoid different ordering between Python 2 and 3.
            mutual_nuclides = natural_nuclides.intersection(library_nuclides)
            absent_nuclides = natural_nuclides.difference(mutual_nuclides)
            mutual_nuclides = sorted(list(mutual_nuclides))
            absent_nuclides = sorted(list(absent_nuclides))

            # If all natural nuclides are present in the library, expand element
            # using all natural nuclides
            if len(absent_nuclides) == 0:
                for nuclide in mutual_nuclides:
                    abundances[nuclide] = NATURAL_ABUNDANCE[nuclide]

            # If no natural elements are present in the library, check if the
            # 0 nuclide is present. If so, set the abundance to 1 for this
            # nuclide. Else, raise an error.
            elif len(mutual_nuclides) == 0:
                nuclide_0 = self.name + '0'
                if nuclide_0 in library_nuclides:
                    abundances[nuclide_0] = 1.0
                else:
                    msg = 'Unable to expand element {0} because the cross '\
                          'section library provided does not contain any of '\
                          'the natural isotopes for that element.'\
                          .format(self.name)
                    raise ValueError(msg)

            # If some, but not all, natural nuclides are in the library, add
            # the mutual nuclides. For the absent nuclides, add them based on
            # our knowledge of the common cross section libraries
            # (ENDF, JEFF, and JENDL)
            else:

                # Add the mutual isotopes
                for nuclide in mutual_nuclides:
                    abundances[nuclide] = NATURAL_ABUNDANCE[nuclide]

                # Adjust the abundances for the absent nuclides
                for nuclide in absent_nuclides:

                    if nuclide in ['O17', 'O18'] and 'O16' in mutual_nuclides:
                        abundances['O16'] += NATURAL_ABUNDANCE[nuclide]
                    elif nuclide == 'Ta180' and 'Ta181' in mutual_nuclides:
                        abundances['Ta181'] += NATURAL_ABUNDANCE[nuclide]
                    elif nuclide == 'W180' and 'W182' in mutual_nuclides:
                        abundances['W182'] += NATURAL_ABUNDANCE[nuclide]
                    else:
                        msg = 'Unsure how to partition natural abundance of ' \
                              'isotope {0} into other natural isotopes of ' \
                              'this element that are present in the cross ' \
                              'section library provided. Consider adding ' \
                              'the isotopes of this element individually.'
                        raise ValueError(msg)

        # If a cross_section library is not present, expand the element into
        # its natural nuclides
        else:
            for nuclide in natural_nuclides:
                abundances[nuclide] = NATURAL_ABUNDANCE[nuclide]

        # Modify mole fractions if enrichment provided
        if enrichment is not None:

            # Calculate the mass fractions of isotopes
            abundances['U234'] = 0.008 * enrichment
            abundances['U235'] = enrichment
            abundances['U238'] = 100.0 - 1.008 * enrichment

            # Convert the mass fractions to mole fractions
            for nuclide in abundances.keys():
                abundances[nuclide] /= atomic_mass(nuclide)

            # Normalize the mole fractions to one
            sum_abundances = sum(abundances.values())
            for nuclide in abundances.keys():
                abundances[nuclide] /= sum_abundances

        # Compute the ratio of the nuclide atomic masses to the element
        # atomic mass
        if percent_type == 'wo':

            # Compute the element atomic mass
            element_am = 0.
            for nuclide in abundances.keys():
                element_am += atomic_mass(nuclide) * abundances[nuclide]

            # Convert the molar fractions to mass fractions
            for nuclide in abundances.keys():
                abundances[nuclide] *= atomic_mass(nuclide) / element_am

            # Normalize the mass fractions to one
            sum_abundances = sum(abundances.values())
            for nuclide in abundances.keys():
                abundances[nuclide] /= sum_abundances

        # Create a list of the isotopes in this element
        isotopes = []
        for nuclide, abundance in abundances.items():
            nuc = openmc.Nuclide(nuclide)
            nuc.scattering = self.scattering
            isotopes.append((nuc, percent * abundance, percent_type))

        return isotopes

    def plot_xs(self, types, divisor_types=None, temperature=294.,
                Erange=(1.E-5, 20.E6), sab_name=None, cross_sections=None,
                enrichment=None, **kwargs):
        """Creates a figure of continuous-energy microscopic cross sections
        for this element

        Parameters
        ----------
        types : Iterable of values of PLOT_TYPES
            The type of cross sections to include in the plot.
        divisor_types : Iterable of values of PLOT_TYPES, optional
            Cross section types which will divide those produced by types
            before plotting. A type of 'unity' can be used to effectively not
            divide some types.
        temperature : float, optional
            Temperature in Kelvin to plot. If not specified, a default
            temperature of 294K will be plotted. Note that the nearest
            temperature in the library for each nuclide will be used as opposed
            to using any interpolation.
        Erange : tuple of floats
            Energy range (in eV) to plot the cross section within
        sab_name : str, optional
            Name of S(a,b) library to apply to MT=2 data when applicable.
        cross_sections : str, optional
            Location of cross_sections.xml file. Default is None.
        enrichment : float, optional
            Enrichment for U235 in weight percent. For example, input 4.95 for
            4.95 weight percent enriched U. Default is None
            (natural composition).
        **kwargs
            All keyword arguments are passed to
            :func:`matplotlib.pyplot.figure`.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Matplotlib Figure of the generated macroscopic cross section

        """

        from matplotlib import pyplot as plt

        E, data = self.calculate_xs(types, temperature, sab_name,
                                    cross_sections)

        if divisor_types:
            cv.check_length('divisor types', divisor_types, len(types),
                            len(types))
            Ediv, data_div = self.calculate_xs(divisor_types, temperature,
                                               sab_name, cross_sections)

            # Create a new union grid, interpolate data and data_div on to that
            # grid, and then do the actual division
            Enum = E[:]
            E = np.union1d(Enum, Ediv)
            data_new = np.zeros((len(types), len(E)))

            for l in range(len(types)):
                data_new[l, :] = \
                    np.divide(np.interp(E, Enum, data[l, :]),
                              np.interp(E, Ediv, data_div[l, :]))
                if divisor_types[l] != 'unity':
                    types[l] = types[l] + ' / ' + divisor_types[l]
            data = data_new

        # Generate the plot
        fig = plt.figure(**kwargs)
        ax = fig.add_subplot(111)
        for i in range(len(data)):
            # Set to loglog or semilogx depending on if we are plotting a data
            # type which we expect to vary linearly
            if types[i] in PLOT_TYPES_LINEAR:
                plot_func = ax.semilogx
            else:
                plot_func = ax.loglog
            if np.sum(data[i, :]) > 0.:
                plot_func(E, data[i, :], label=types[i])

        ax.set_xlabel('Energy [eV]')
        if divisor_types:
            ax.set_ylabel('Elemental Cross Section')
        else:
            ax.set_ylabel('Elemental Cross Section [b]')
        ax.legend(loc='best')
        ax.set_xlim(Erange)
        if self.name is not None:
            title = 'Cross Section for ' + self.name
            ax.set_title(title)

        return fig

    def calculate_xs(self, types, temperature=294., sab_name=None,
                     cross_sections=None, enrichment=None):
        """Calculates continuous-energy macroscopic cross sections of a
        requested type

        Parameters
        ----------
        types : Iterable of values of PLOT_TYPES
            The type of cross sections to calculate
        temperature : float, optional
            Temperature in Kelvin to plot. If not specified, a default
            temperature of 294K will be plotted. Note that the nearest
            temperature in the library for each nuclide will be used as opposed
            to using any interpolation.
        sab_name : str, optional
            Name of S(a,b) library to apply to MT=2 data when applicable.
        cross_sections : str, optional
            Location of cross_sections.xml file. Default is None.
        enrichment : float, optional
            Enrichment for U235 in weight percent. For example, input 4.95 for
            4.95 weight percent enriched U. Default is None
            (natural composition).

        Returns
        -------
        unionE : numpy.array
            Energies at which cross sections are calculated, in units of eV
        data : numpy.ndarray
            Macroscopic cross sections calculated at the energy grid described
            by unionE

        """

        # Check types
        if cross_sections is not None:
            cv.check_type('cross_sections', cross_sections, str)
        cv.check_iterable_type('types', types, str)

        # Parse the types
        mts = []
        ops = []
        yields = []
        for line in types:
            if line in PLOT_TYPES:
                mts.append(PLOT_TYPES_MT[line])
                yields.append(PLOT_TYPES_YIELD[line])
                ops.append(PLOT_TYPES_OP[line])
            else:
                # Not a built-in type, we have to parse it ourselves
                raise NotImplementedError()

        cv.check_type('temperature', temperature, Real)

        # If cross_sections is None, get the cross sections from the
        # OPENMC_CROSS_SECTIONS environment variable
        if cross_sections is None:
            cross_sections = os.environ.get('OPENMC_CROSS_SECTIONS')

        # If a cross_sections library is present, check natural nuclides
        # against the nuclides in the library
        if cross_sections is not None:
            library = openmc.data.DataLibrary.from_xml(cross_sections)
        else:
            raise ValueError("cross_sections or OPENMC_CROSS_SECTIONS "
                             "environmental variable must be set")

        # Expand elements in to nuclides with atomic densities
        nuclides = self.expand(100., 'ao', enrichment=enrichment,
                               cross_sections=cross_sections)

        # For ease of processing split out nuc and nuc_density
        nuc_fractions = [nuclide[1] for nuclide in nuclides]

        # Identify the nuclides which have S(a,b) data
        sabs = {}
        for nuclide in nuclides:
            sabs[nuclide[0].name] = None
        if sab_name:
            sab = openmc.data.ThermalScattering.from_hdf5(sab_name)
            for nuc in sab.nuclides:
                sabs[nuc] = library.get_by_material(sab_name)['path']

        # Now we can create the data sets to be plotted
        xs = []
        E = []
        for nuclide in nuclides:
            sab_tab = sabs[nuclide[0].name]
            temp_E, temp_xs = nuclide[0].calculate_xs(types, temperature,
                                                      sab_tab, cross_sections)
            E.append(temp_E)
            xs.append(temp_xs)

        # Condense the data for every nuclide
        # First create a union energy grid
        unionE = E[0]
        for n in range(1, len(E)):
            unionE = np.union1d(unionE, E[n])

        # Now we can combine all the nuclidic data
        data = np.zeros((len(mts), len(unionE)))
        for l in range(len(mts)):
            if types[l] == 'unity':
                data[l, :] = 1.
            else:
                for n in range(len(nuclides)):
                    data[l, :] += nuc_fractions[n] * xs[n][l](unionE)

        return unionE, data
