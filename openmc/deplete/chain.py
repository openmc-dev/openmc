"""chain module.

This module contains information about a depletion chain.  A depletion chain is
loaded from an .xml file and all the nuclides are linked together.
"""

from collections import OrderedDict, defaultdict
from io import StringIO
from itertools import chain
import math
import re
import os

# Try to use lxml if it is available. It preserves the order of attributes and
# provides a pretty-printer by default. If not available, use OpenMC function to
# pretty print.
try:
    import lxml.etree as ET
    _have_lxml = True
except ImportError:
    import xml.etree.ElementTree as ET
    _have_lxml = False
import scipy.sparse as sp

import openmc.data
from openmc.clean_xml import clean_xml_indentation
from .nuclide import Nuclide, DecayTuple, ReactionTuple


# tuple of (reaction name, possible MT values, (dA, dZ)) where dA is the change
# in the mass number and dZ is the change in the atomic number
_REACTIONS = [
    ('(n,2n)', set(chain([16], range(875, 892))), (-1, 0)),
    ('(n,3n)', {17}, (-2, 0)),
    ('(n,4n)', {37}, (-3, 0)),
    ('(n,gamma)', {102}, (1, 0)),
    ('(n,p)', set(chain([103], range(600, 650))), (0, -1)),
    ('(n,a)', set(chain([107], range(800, 850))), (-3, -2))
]


def replace_missing(product, decay_data):
    """Replace missing product with suitable decay daughter.

    Parameters
    ----------
    product : str
        Name of product in GND format, e.g. 'Y86_m1'.
    decay_data : dict
        Dictionary of decay data

    Returns
    -------
    product : str
        Replacement for missing product in GND format.

    """
    # Determine atomic number, mass number, and metastable state
    Z, A, state = openmc.data.zam(product)
    symbol = openmc.data.ATOMIC_SYMBOL[Z]

    # Replace neutron with proton
    if Z == 0 and A == 1:
        return 'H1'

    # First check if ground state is available
    if state:
        product = '{}{}'.format(symbol, A)

    # Find isotope with longest half-life
    half_life = 0.0
    for nuclide, data in decay_data.items():
        m = re.match(r'{}(\d+)(?:_m\d+)?'.format(symbol), nuclide)
        if m:
            # If we find a stable nuclide, stop search
            if data.nuclide['stable']:
                mass_longest_lived = int(m.group(1))
                break
            if data.half_life.nominal_value > half_life:
                mass_longest_lived = int(m.group(1))
                half_life = data.half_life.nominal_value

    # If mass number of longest-lived isotope is less than that of missing
    # product, assume it undergoes beta-. Otherwise assume beta+.
    beta_minus = (mass_longest_lived < A)

    # Iterate until we find an existing nuclide
    while product not in decay_data:
        if Z > 98:
            Z -= 2
            A -= 4
        else:
            if beta_minus:
                Z += 1
            else:
                Z -= 1
        product = '{}{}'.format(openmc.data.ATOMIC_SYMBOL[Z], A)

    return product


class Chain(object):
    """Full representation of a depletion chain.

    A depletion chain can be created by using the :meth:`from_endf` method which
    requires a list of ENDF incident neutron, decay, and neutron fission product
    yield sublibrary files. The depletion chain used during a depletion
    simulation is indicated by either an argument to
    :class:`openmc.deplete.Operator` or through the
    :envvar:`OPENMC_DEPLETE_CHAIN` environment variable.

    Attributes
    ----------
    nuclides : list of openmc.deplete.Nuclide
        Nuclides present in the chain.
    reactions : list of str
        Reactions that are tracked in the depletion chain
    nuclide_dict : OrderedDict of str to int
        Maps a nuclide name to an index in nuclides.

    """

    def __init__(self):
        self.nuclides = []
        self.reactions = []
        self.nuclide_dict = OrderedDict()

    def __contains__(self, nuclide):
        return nuclide in self.nuclide_dict

    def __getitem__(self, name):
        """Get a Nuclide by name."""
        return self.nuclides[self.nuclide_dict[name]]

    def __len__(self):
        """Number of nuclides in chain."""
        return len(self.nuclides)

    @classmethod
    def from_endf(cls, decay_files, fpy_files, neutron_files):
        """Create a depletion chain from ENDF files.

        Parameters
        ----------
        decay_files : list of str
            List of ENDF decay sub-library files
        fpy_files : list of str
            List of ENDF neutron-induced fission product yield sub-library files
        neutron_files : list of str
            List of ENDF neutron reaction sub-library files

        """
        chain = cls()

        # Create dictionary mapping target to filename
        print('Processing neutron sub-library files...')
        reactions = {}
        for f in neutron_files:
            evaluation = openmc.data.endf.Evaluation(f)
            name = evaluation.gnd_name
            reactions[name] = {}
            for mf, mt, nc, mod in evaluation.reaction_list:
                if mf == 3:
                    file_obj = StringIO(evaluation.section[3, mt])
                    openmc.data.endf.get_head_record(file_obj)
                    q_value = openmc.data.endf.get_cont_record(file_obj)[1]
                    reactions[name][mt] = q_value

        # Determine what decay and FPY nuclides are available
        print('Processing decay sub-library files...')
        decay_data = {}
        for f in decay_files:
            data = openmc.data.Decay(f)
            # Skip decay data for neutron itself
            if data.nuclide['atomic_number'] == 0:
                continue
            decay_data[data.nuclide['name']] = data

        print('Processing fission product yield sub-library files...')
        fpy_data = {}
        for f in fpy_files:
            data = openmc.data.FissionProductYields(f)
            fpy_data[data.nuclide['name']] = data

        print('Creating depletion_chain...')
        missing_daughter = []
        missing_rx_product = []
        missing_fpy = []
        missing_fp = []

        for idx, parent in enumerate(sorted(decay_data, key=openmc.data.zam)):
            data = decay_data[parent]

            nuclide = Nuclide()
            nuclide.name = parent

            chain.nuclides.append(nuclide)
            chain.nuclide_dict[parent] = idx

            if not data.nuclide['stable'] and data.half_life.nominal_value != 0.0:
                nuclide.half_life = data.half_life.nominal_value
                nuclide.decay_energy = sum(E.nominal_value for E in
                                           data.average_energies.values())
                sum_br = 0.0
                for i, mode in enumerate(data.modes):
                    type_ = ','.join(mode.modes)
                    if mode.daughter in decay_data:
                        target = mode.daughter
                    else:
                        print('missing {} {} {}'.format(parent, ','.join(mode.modes), mode.daughter))
                        target = replace_missing(mode.daughter, decay_data)

                    # Write branching ratio, taking care to ensure sum is unity
                    br = mode.branching_ratio.nominal_value
                    sum_br += br
                    if i == len(data.modes) - 1 and sum_br != 1.0:
                        br = 1.0 - sum(m.branching_ratio.nominal_value
                                       for m in data.modes[:-1])

                    # Append decay mode
                    nuclide.decay_modes.append(DecayTuple(type_, target, br))

            if parent in reactions:
                reactions_available = set(reactions[parent].keys())
                for name, mts, changes in _REACTIONS:
                    if mts & reactions_available:
                        delta_A, delta_Z = changes
                        A = data.nuclide['mass_number'] + delta_A
                        Z = data.nuclide['atomic_number'] + delta_Z
                        daughter = '{}{}'.format(openmc.data.ATOMIC_SYMBOL[Z], A)

                        if name not in chain.reactions:
                            chain.reactions.append(name)

                        if daughter not in decay_data:
                            missing_rx_product.append((parent, name, daughter))

                        # Store Q value
                        for mt in sorted(mts):
                            if mt in reactions[parent]:
                                q_value = reactions[parent][mt]
                                break
                        else:
                            q_value = 0.0

                        nuclide.reactions.append(ReactionTuple(
                            name, daughter, q_value, 1.0))

                if any(mt in reactions_available for mt in [18, 19, 20, 21, 38]):
                    if parent in fpy_data:
                        q_value = reactions[parent][18]
                        nuclide.reactions.append(
                            ReactionTuple('fission', 0, q_value, 1.0))

                        if 'fission' not in chain.reactions:
                            chain.reactions.append('fission')
                    else:
                        missing_fpy.append(parent)

            if parent in fpy_data:
                fpy = fpy_data[parent]

                if fpy.energies is not None:
                    nuclide.yield_energies = fpy.energies
                else:
                    nuclide.yield_energies = [0.0]

                for E, table in zip(nuclide.yield_energies, fpy.independent):
                    yield_replace = 0.0
                    yields = defaultdict(float)
                    for product, y in table.items():
                        # Handle fission products that have no decay data available
                        if product not in decay_data:
                            daughter = replace_missing(product, decay_data)
                            product = daughter
                            yield_replace += y.nominal_value

                        yields[product] += y.nominal_value

                    if yield_replace > 0.0:
                        missing_fp.append((parent, E, yield_replace))

                    nuclide.yield_data[E] = []
                    for k in sorted(yields, key=openmc.data.zam):
                        nuclide.yield_data[E].append((k, yields[k]))

        # Display warnings
        if missing_daughter:
            print('The following decay modes have daughters with no decay data:')
            for mode in missing_daughter:
                print('  {}'.format(mode))
            print('')

        if missing_rx_product:
            print('The following reaction products have no decay data:')
            for vals in missing_rx_product:
                print('{} {} -> {}'.format(*vals))
            print('')

        if missing_fpy:
            print('The following fissionable nuclides have no fission product yields:')
            for parent in missing_fpy:
                print('  ' + parent)
            print('')

        if missing_fp:
            print('The following nuclides have fission products with no decay data:')
            for vals in missing_fp:
                print('  {}, E={} eV (total yield={})'.format(*vals))

        return chain

    @classmethod
    def from_xml(cls, filename):
        """Reads a depletion chain XML file.

        Parameters
        ----------
        filename : str
            The path to the depletion chain XML file.

        """
        chain = cls()

        # Load XML tree
        root = ET.parse(str(filename))

        for i, nuclide_elem in enumerate(root.findall('nuclide')):
            nuc = Nuclide.from_xml(nuclide_elem)
            chain.nuclide_dict[nuc.name] = i

            # Check for reaction paths
            for rx in nuc.reactions:
                if rx.type not in chain.reactions:
                    chain.reactions.append(rx.type)

            chain.nuclides.append(nuc)

        return chain

    def export_to_xml(self, filename):
        """Writes a depletion chain XML file.

        Parameters
        ----------
        filename : str
            The path to the depletion chain XML file.

        """

        root_elem = ET.Element('depletion_chain')
        for nuclide in self.nuclides:
            root_elem.append(nuclide.to_xml_element())

        tree = ET.ElementTree(root_elem)
        if _have_lxml:
            tree.write(str(filename), encoding='utf-8', pretty_print=True)
        else:
            clean_xml_indentation(root_elem)
            tree.write(str(filename), encoding='utf-8')

    def form_matrix(self, rates):
        """Forms depletion matrix.

        Parameters
        ----------
        rates : numpy.ndarray
            2D array indexed by (nuclide, reaction)

        Returns
        -------
        scipy.sparse.csr_matrix
            Sparse matrix representing depletion.

        """
        matrix = defaultdict(float)
        reactions = set()

        for i, nuc in enumerate(self.nuclides):

            if nuc.n_decay_modes != 0:
                # Decay paths
                # Loss
                decay_constant = math.log(2) / nuc.half_life

                if decay_constant != 0.0:
                    matrix[i, i] -= decay_constant

                # Gain
                for _, target, branching_ratio in nuc.decay_modes:
                    # Allow for total annihilation for debug purposes
                    if target != 'Nothing':
                        branch_val = branching_ratio * decay_constant

                        if branch_val != 0.0:
                            k = self.nuclide_dict[target]
                            matrix[k, i] += branch_val

            if nuc.name in rates.index_nuc:
                # Extract all reactions for this nuclide in this cell
                nuc_ind = rates.index_nuc[nuc.name]
                nuc_rates = rates[nuc_ind, :]

                for r_type, target, _, br in nuc.reactions:
                    # Extract reaction index, and then final reaction rate
                    r_id = rates.index_rx[r_type]
                    path_rate = nuc_rates[r_id]

                    # Loss term -- make sure we only count loss once for
                    # reactions with branching ratios
                    if r_type not in reactions:
                        reactions.add(r_type)
                        if path_rate != 0.0:
                            matrix[i, i] -= path_rate

                    # Gain term; allow for total annihilation for debug purposes
                    if target != 'Nothing':
                        if r_type != 'fission':
                            if path_rate != 0.0:
                                k = self.nuclide_dict[target]
                                matrix[k, i] += path_rate * br
                        else:
                            # Assume that we should always use thermal fission
                            # yields. At some point it would be nice to account
                            # for the energy-dependence..
                            energy, data = sorted(nuc.yield_data.items())[0]
                            for product, y in data:
                                yield_val = y * path_rate
                                if yield_val != 0.0:
                                    k = self.nuclide_dict[product]
                                    matrix[k, i] += yield_val

                # Clear set of reactions
                reactions.clear()

        # Use DOK matrix as intermediate representation, then convert to CSR and return
        n = len(self)
        matrix_dok = sp.dok_matrix((n, n))
        dict.update(matrix_dok, matrix)
        return matrix_dok.tocsr()
