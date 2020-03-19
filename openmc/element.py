from collections import OrderedDict
import re
import os
from xml.etree import ElementTree as ET

import openmc.checkvalue as cv
from numbers import Real
from openmc.data import NATURAL_ABUNDANCE, atomic_mass


class Element(str):
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

    """

    def __new__(cls, name):
        cv.check_type('element name', name, str)
        cv.check_length('element name', name, 1, 2)
        return super().__new__(cls, name)

    @property
    def name(self):
        return self

    def expand(self, percent, percent_type, enrichment=None,
               enrichment_target=None, enrichment_type=None,
               cross_sections=None):
        """Expand natural element into its naturally-occurring isotopes.

        An optional cross_sections argument or the :envvar:`OPENMC_CROSS_SECTIONS`
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
            Enrichment of an enrichment_taget nuclide in percent (ao or wo).
            If enrichment_taget is not supplied then it is enrichment for U235
            in weight percent. For example, input 4.95 for 4.95 weight percent
            enriched U. Default is None (natural composition).
        enrichment_target: str, optional
            Single nuclide name to enrich from a natural composition (e.g., 'O16')
        enrichment_type: {'ao', 'wo'}, optional
            'ao' for enrichment as atom percent and 'wo' for weight percent.
            Default is: 'ao' for two-isotope enrichment; 'wo' for U enrichment
        cross_sections : str, optional
            Location of cross_sections.xml file. Default is None.

        Returns
        -------
        isotopes : list
            Naturally-occurring isotopes of the element. Each item of the list
            is a tuple consisting of a nuclide string, the atom/weight percent,
            and the string 'ao' or 'wo'.

        Raises
        ------
        ValueError
            No data is available for any of natural isotopes of the element
        ValueError
            If only some natural isotopes are available in the cross-section data
            library and the element is not O, W, or Ta
        ValueError
            If a non-naturally-occurring isotope is requested
        ValueError
            If enrichment is requested of an element with more than two
            naturally-occurring isotopes.
        ValueError
            If enrichment procedure for Uranium is used when element is not
            Uranium.
        ValueError
            Uranium enrichment is requested with enrichment_type=='ao'

        Notes
        -----
        When the `enrichment` argument is specified, a correlation from
        `ORNL/CSD/TM-244 <https://doi.org/10.2172/5561567>`_ is used to
        calculate the weight fractions of U234, U235, U236, and U238. Namely,
        the weight fraction of U234 and U236 are taken to be 0.89% and 0.46%,
        respectively, of the U235 weight fraction. The remainder of the
        isotopic weight is assigned to U238.

        When the `enrichment` argument is specified with `enrichment_target`, a
        general enrichment procedure is used for elements composed of exactly
        two naturally-occurring isotopes. `enrichment` is interpreted as atom
        percent by default but can be controlled by the `enrichment_type`
        argument.

        """
        # Check input
        if enrichment_type is not None:
            cv.check_value('enrichment_type', enrichment_type, {'ao', 'wo'})

        if enrichment is not None:
            cv.check_less_than('enrichment', enrichment, 100.0, equality=True)
            cv.check_greater_than('enrichment', enrichment, 0., equality=True)

        # Get the nuclides present in nature
        natural_nuclides = set()
        for nuclide in sorted(NATURAL_ABUNDANCE.keys()):
            if re.match(r'{}\d+'.format(self), nuclide):
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
            for child in root.findall('library'):
                nuclide = child.attrib['materials']
                if re.match(r'{}\d+'.format(self), nuclide) and \
                   '_m' not in nuclide:
                    library_nuclides.add(nuclide)

            # Get a set of the mutual and absent nuclides. Convert to lists
            # and sort to avoid different ordering between Python 2 and 3.
            mutual_nuclides = natural_nuclides.intersection(library_nuclides)
            absent_nuclides = natural_nuclides.difference(mutual_nuclides)
            mutual_nuclides = sorted(list(mutual_nuclides))
            absent_nuclides = sorted(list(absent_nuclides))

            # If all natural nuclides are present in the library,
            # expand element using all natural nuclides
            if len(absent_nuclides) == 0:
                for nuclide in mutual_nuclides:
                    abundances[nuclide] = NATURAL_ABUNDANCE[nuclide]

            # If no natural elements are present in the library, check if the
            # 0 nuclide is present. If so, set the abundance to 1 for this
            # nuclide. Else, raise an error.
            elif len(mutual_nuclides) == 0:
                nuclide_0 = self + '0'
                if nuclide_0 in library_nuclides:
                    abundances[nuclide_0] = 1.0
                else:
                    msg = 'Unable to expand element {0} because the cross '\
                          'section library provided does not contain any of '\
                          'the natural isotopes for that element.'\
                          .format(self)
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
        # Old treatment for Uranium
        if enrichment is not None and enrichment_target is None:

            # Check that the element is Uranium
            if self.name != 'U':
                msg = ('Enrichment procedure for Uranium was requested, '
                       'but the isotope is {} not U'.format(self))
                raise ValueError(msg)

            # Check that enrichment_type is not 'ao'
            if enrichment_type == 'ao':
                msg = ('Enrichment procedure for Uranium requires that '
                       'enrichment value is provided as wo%.')
                raise ValueError(msg)

            # Calculate the mass fractions of isotopes
            abundances['U234'] = 0.0089 * enrichment
            abundances['U235'] = enrichment
            abundances['U236'] = 0.0046 * enrichment
            abundances['U238'] = 100.0 - 1.0135 * enrichment

            # Convert the mass fractions to mole fractions
            for nuclide in abundances.keys():
                abundances[nuclide] /= atomic_mass(nuclide)

            # Normalize the mole fractions to one
            sum_abundances = sum(abundances.values())
            for nuclide in abundances.keys():
                abundances[nuclide] /= sum_abundances

        # Modify mole fractions if enrichment provided
        # New treatment for arbitrary element
        elif enrichment is not None and enrichment_target is not None:

            # Provide more informative error message for U235
            if enrichment_target == 'U235':
                msg = ("There is a special procedure for enrichment of U235 "
                       "in U. To invoke it, the arguments 'enrichment_target'"
                       "and 'enrichment_type' should be omitted. Provide "
                       "a value only for 'enrichment' in weight percent.")
                raise ValueError(msg)

            # Check if it is two-isotope mixture
            if len(abundances) != 2:
                msg = ('Element {} does not consist of two naturally-occurring '
                       'isotopes. Please enter isotopic abundances manually.'
                       .format(self))
                raise ValueError(msg)

            # Check if the target nuclide is present in the mixture
            if enrichment_target not in abundances:
                msg = ('The target nuclide {} is not one of the naturally-occurring '
                       'isotopes ({})'.format(enrichment_target, list(abundances)))
                raise ValueError(msg)

            # If weight percent enrichment is requested convert to mass fractions
            if enrichment_type == 'wo':
                # Convert the atomic abundances to weight fractions
                # Compute the element atomic mass
                element_am = sum(atomic_mass(nuc)*abundances[nuc] for nuc in abundances)

                # Convert Molar Fractions to mass fractions
                for nuclide in abundances:
                    abundances[nuclide] *= atomic_mass(nuclide) / element_am

                # Normalize to one
                sum_abundances = sum(abundances.values())
                for nuclide in abundances:
                    abundances[nuclide] /= sum_abundances

            # Enrich the mixture
            # The procedure is more generic that it needs to be. It allows
            # to enrich mixtures of more then 2 isotopes, keeping the ratios
            # of non-enriched nuclides the same as in natural composition

            # Get fraction of non-enriched isotopes in nat. composition
            non_enriched = 1.0 - abundances[enrichment_target]
            tail_fraction = 1.0 - enrichment / 100.0

            # Enrich all nuclides
            # Do bogus operation for enrichment target but overwrite immediatly
            # to avoid if statement in the loop
            for nuclide, fraction in abundances.items():
                abundances[nuclide] = tail_fraction * fraction / non_enriched
            abundances[enrichment_target] = enrichment / 100.0

            # Convert back to atomic fractions if requested
            if enrichment_type == 'wo':
                # Convert the mass fractions to mole fractions
                for nuclide in abundances:
                    abundances[nuclide] /= atomic_mass(nuclide)

                # Normalize the mole fractions to one
                sum_abundances = sum(abundances.values())
                for nuclide in abundances:
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
            isotopes.append((nuclide, percent * abundance, percent_type))

        return isotopes
