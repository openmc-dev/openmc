import re
import sys
import os

from six import string_types

import openmc
from openmc.checkvalue import check_type, check_length
from openmc.data import NATURAL_ABUNDANCE
from xml.etree import ElementTree as ET

class Element(object):
    """A natural element used in a material via <element>. Internally, OpenMC
    will expand the natural element into isotopes based on the known natural
    abundances.

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
        check_type('element name', name, string_types)
        check_length('element name', name, 1, 2)
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
            Enrichment percent for U235 in U. Default is None
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

        # Get the length of this elements atomic symbol
        name_len = len(self.name)

        # Get the nuclides present in nature
        natural_nuclides = set()
        for nuclide in sorted(NATURAL_ABUNDANCE.keys()):
            if re.match(r'{}\d+'.format(self.name), nuclide):
                natural_nuclides.add(int(nuclide[name_len:]))

        # Create lists to store the expanded nuclides and abundances
        nuclides = []
        abundances = []

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
                    library_nuclides.add(int(nuclide[name_len:]))

            # Get a set of the mutual and absent nuclides
            mutual_nuclides = natural_nuclides.intersection(library_nuclides)
            absent_nuclides = natural_nuclides.difference(mutual_nuclides)
            mutual_nuclides = sorted(list(mutual_nuclides))
            absent_nuclides = sorted(list(absent_nuclides))

            # If all natural nuclides are present in the library, expand element
            # using all natural nuclides
            if len(absent_nuclides) == 0:
                for nuclide, abundance in sorted(NATURAL_ABUNDANCE.items()):
                    if re.match(r'{}\d+'.format(self.name), nuclide):
                        nuc = openmc.Nuclide(nuclide)
                        nuc.scattering = self.scattering
                        nuclides.append(nuc)
                        abundances.append(abundance)

            # If no natural elements are present in the library, check if the
            # 0 element is present. If so, set the abundance to 1 for this
            # nuclide. Else, raise an error.
            elif len(mutual_nuclides) == 0:
                if 0 in library_nuclides:
                    nuc = openmc.Nuclide(self.name + '0')
                    nuc.scattering = self.scattering
                    nuclides.append(nuc)
                    abundances.append(1.0)
                else:
                    msg = 'Unable to expand element {0} because the cross '\
                          'section library provided does not contain any of '\
                          'the natural isotopes for that element.'\
                          .format(self.name)
                    raise ValueError(msg)

            # If some, but not all, natural nuclides are in the library, add
            # the mutual nuclides. For the absent nuclides, increment the
            # abundance of the nearest mutual nuclide with their abundances.
            else:

                # Add the mutual isotopes
                nuclides_a = []
                for nuclide, abundance in sorted(NATURAL_ABUNDANCE.items()):
                    if re.match(r'{}\d+'.format(self.name), nuclide) and \
                       int(nuclide[name_len:]) in mutual_nuclides:
                        nuc = openmc.Nuclide(nuclide)
                        nuc.scattering = self.scattering
                        nuclides.append(nuc)
                        nuclides_a.append(int(nuclide[name_len:]))
                        abundances.append(abundance)

                # Adjust the abundances for the absent nuclides
                for nuclide, abundance in sorted(NATURAL_ABUNDANCE.items()):
                    if re.match(r'{}\d+'.format(self.name), nuclide) and \
                       int(nuclide[name_len:]) in absent_nuclides:

                        # Get index to the nearest nuclide
                        a = int(nuclide[name_len:])
                        i = min(list(range(len(nuclides_a))), key=lambda j: \
                                abs(nuclides_a[j] - a))

                        # Increment abundance of the nearest nuclide
                        abundances[i] += abundance

        # If a cross_section library is not present, expand the element into
        # its natural nuclides
        else:
            for nuclide, abundance in sorted(NATURAL_ABUNDANCE.items()):
                if re.match(r'{}\d+'.format(self.name), nuclide):
                    nuclides.append(openmc.Nuclide(nuclide))
                    abundances.append(abundance)

        # Create a list of atomic masses
        n_nuclides = len(nuclides)
        atomic_masses = []
        for nuclide in nuclides:
            atomic_masses.append(openmc.data.atomic_mass(nuclide.name))

        # Modify mole fractions if enrichment provided
        if enrichment is not None:

            # Get the indices for the uranium nuclides
            for i,nuc in enumerate(nuclides):
                if nuc.name == 'U234':
                    u234 = i
                elif nuc.name == 'U235':
                    u235 = i
                elif nuc.name == 'U238':
                    u238 = i

            # Calculate the mass fractions of isotopes
            abundances[u234] = 0.008 * enrichment
            abundances[u235] = enrichment
            abundances[u238] = 1.0 - 1.008 * enrichment

            # Convert the mass fractions to mole fractions
            for i in range(n_nuclides):
                abundances[i] /= atomic_masses[i]

            # Normalize the mole fractions to one
            sum_abundances = sum(abundances)
            for i in range(n_nuclides):
                abundances[i] /= sum_abundances

        # Compute the ratio of the nuclide atomic massess to the element
        # atomic mass
        if percent_type == 'wo':

            # Compute the element atomic mass
            element_am = 0.
            for i in range(n_nuclides):
                element_am += atomic_masses[i] * abundances[i]

            # Convert the molar fractions to mass fractions
            for i in range(n_nuclides):
                abundances[i] *= atomic_masses[i] / element_am

            # Normalize the mass fractions to one
            sum_abundances = sum(abundances)
            for i in range(n_nuclides):
                abundances[i] /= sum_abundances

        # Create a list of the isotopes in this element
        isotopes = []
        for nuclide, abundance in zip(nuclides, abundances):
            #pct = float('{:2.8g}'.format())
            pct = percent*abundance
            isotopes.append((nuclide, pct, percent_type))

        return isotopes
