"""function module.

This module contains the Operator class, which is then passed to an integrator
to run a full depletion simulation.
"""

from collections import namedtuple
import os
from pathlib import Path
from abc import ABC, abstractmethod
from xml.etree import ElementTree as ET
from warnings import warn
from numbers import Real

from numpy import nonzero, empty

from openmc.data import DataLibrary, JOULE_PER_EV
from openmc.checkvalue import check_type, check_greater_than
from .chain import Chain

OperatorResult = namedtuple('OperatorResult', ['k', 'rates'])
OperatorResult.__doc__ = """\
Result of applying transport operator

Parameters
----------
k : float
    Resulting eigenvalue
rates : openmc.deplete.ReactionRates
    Resulting reaction rates

"""
try:
    OperatorResult.k.__doc__ = None
    OperatorResult.rates.__doc__ = None
except AttributeError:
    # Can't set __doc__ on properties on Python 3.4
    pass


class TransportOperator(ABC):
    """Abstract class defining a transport operator

    Each depletion integrator is written to work with a generic transport
    operator that takes a vector of material compositions and returns an
    eigenvalue and reaction rates. This abstract class sets the requirements for
    such a transport operator. Users should instantiate
    :class:`openmc.deplete.Operator` rather than this class.

    Parameters
    ----------
    chain_file : str, optional
        Path to the depletion chain XML file.  Defaults to the file
        listed under ``depletion_chain`` in
        :envvar:`OPENMC_CROSS_SECTIONS` environment variable.
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``.
    dilute_initial : float, optional
        Initial atom density [atoms/cm^3] to add for nuclides that are zero
        in initial condition to ensure they exist in the decay chain.
        Only done for nuclides with reaction rates.
        Defaults to 1.0e3.

    Attributes
    ----------
    dilute_initial : float
        Initial atom density [atoms/cm^3] to add for nuclides that are zero
        in initial condition to ensure they exist in the decay chain.
        Only done for nuclides with reaction rates.
    """
    def __init__(self, chain_file=None, fission_q=None, dilute_initial=1.0e3):
        self.dilute_initial = dilute_initial
        self.output_dir = '.'

        # Read depletion chain
        if chain_file is None:
            chain_file = os.environ.get("OPENMC_DEPLETE_CHAIN", None)
            if chain_file is None:
                data = DataLibrary.from_xml()
                # search for depletion_chain path from end of list
                for lib in reversed(data.libraries):
                    if lib['type'] == 'depletion_chain':
                        break
                else:
                    raise IOError(
                        "No chain specified, either manually or "
                        "under depletion_chain in environment variable "
                        "OPENMC_CROSS_SECTIONS.")
                chain_file = lib['path']
            else:
                warn("Use of OPENMC_DEPLETE_CHAIN is deprecated in favor "
                     "of adding depletion_chain to OPENMC_CROSS_SECTIONS",
                     FutureWarning)
        self.chain = Chain.from_xml(chain_file, fission_q)

    @property
    def dilute_initial(self):
        """Initial atom density for nuclides with zero initial concentration"""
        return self._dilute_initial

    @dilute_initial.setter
    def dilute_initial(self, value):
        check_type("dilute_initial", value, Real)
        check_greater_than("dilute_initial", value, 0.0, equality=True)
        self._dilute_initial = value

    @abstractmethod
    def __call__(self, vec, print_out=True):
        """Runs a simulation.

        Parameters
        ----------
        vec : list of numpy.ndarray
            Total atoms to be used in function.
        print_out : bool, optional
            Whether or not to print out time.

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """
        pass

    def __enter__(self):
        # Save current directory and move to specific output directory
        self._orig_dir = os.getcwd()
        if not self.output_dir.exists():
            self.output_dir.mkdir()  # exist_ok parameter is 3.5+

        # In Python 3.6+, chdir accepts a Path directly
        os.chdir(str(self.output_dir))

        return self.initial_condition()

    def __exit__(self, exc_type, exc_value, traceback):
        self.finalize()
        os.chdir(self._orig_dir)

    @property
    def output_dir(self):
        return self._output_dir

    @output_dir.setter
    def output_dir(self, output_dir):
        self._output_dir = Path(output_dir)

    @abstractmethod
    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.ndarray
            Total density for initial conditions.
        """

        pass

    @abstractmethod
    def get_results_info(self):
        """Returns volume list, cell lists, and nuc lists.

        Returns
        -------
        volume : list of float
            Volumes corresponding to materials in burn_list
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all cell IDs to be burned.  Used for sorting the simulation.
        full_burn_list : list of int
            All burnable materials in the geometry.
        """

        pass

    def finalize(self):
        pass


class ReactionRateHelper(ABC):
    """Abstract class for generating reaction rates for operators

    Responsible for generating reaction rate tallies for burnable
    materials, given nuclides and scores from the operator.

    Reaction rates are passed back to the operator for be used in
    an :class:`openmc.deplete.OperatorResult` instance

    Parameters
    ----------
    n_nucs : int
        Number of burnable nuclides tracked by :class:`openmc.deplete.Operator`
    n_react : int
        Number of reactions tracked by :class:`openmc.deplete.Operator`

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates.
    """

    def __init__(self, n_nucs, n_react):
        self._nuclides = None
        self._rate_tally = None
        self._results_cache = empty((n_nucs, n_react))

    @abstractmethod
    def generate_tallies(self, materials, scores):
        """Use the C API to build tallies needed for reaction rates"""

    @property
    def nuclides(self):
        """List of nuclides with requested reaction rates"""
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        check_type("nuclides", nuclides, list, str)
        self._nuclides = nuclides
        self._rate_tally.nuclides = nuclides

    @abstractmethod
    def get_material_rates(self, mat_id, nuc_index, react_index):
        """Return 2D array of [nuclide, reaction] reaction rates

        Parameters
        ----------
        mat_id : int
            Unique ID for the requested material
        nuc_index : list of str
            Ordering of desired nuclides
        react_index : list of str
            Ordering of reactions
        """

    def divide_by_adens(self, number):
        """Normalize reaction rates by number of nuclides

        Acts on the current material examined by
        :meth:`get_material_rates`

        Parameters
        ----------
        number : iterable of float
            Number density [atoms/b-cm] of each nuclide tracked in the calculation.

        Returns
        -------
        results : numpy.ndarray
            Array of reactions rates of shape ``(n_nuclides, n_rxns)``
            normalized by the number of nuclides
        """

        mask = nonzero(number)
        results = self._results_cache
        for col in range(results.shape[1]):
            results[mask, col] /= number[mask]
        return results


class EnergyHelper(ABC):
    """Abstract class for obtaining energy produced

    The ultimate goal of this helper is to provide instances of
    :class:`openmc.deplete.Operator` with the total energy produced
    in a transport simulation. This information, provided with the
    power requested by the user and reaction rates from a
    :class:`ReactionRateHelper` will scale reaction rates to the
    correct values.

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates. Ordered to be
        consistent with :class:`openmc.deplete.Operator`
    energy : float
        Total energy [J/s/source neutron] produced in a transport simulation.
        Updated in the material iteration with :meth:`update`.
    """

    def __init__(self):
        self._nuclides = None
        self._energy = 0.0

    @property
    def energy(self):
        return self._energy * JOULE_PER_EV

    def reset(self):
        """Reset energy produced prior to unpacking tallies"""
        self._energy = 0.0

    @abstractmethod
    def prepare(self, chain_nucs, rate_index, materials):
        """Perform work needed to obtain energy produced

        This method is called prior to the transport simulations
        in :meth:`openmc.deplete.Operator.initial_condition`.

        Parameters
        ----------
        chain_nucs : list of str
            All nuclides to be tracked in this problem
        rate_index : dict of str to int
            Mapping from nuclide name to index in the
            `fission_rates` for :meth:`update`.
        materials : list of str
            All materials tracked on the operator helped by this
            object. Should correspond to
            :attr:`openmc.deplete.Operator.burnable_materials`
        """

    def update(self, fission_rates, mat_index):
        """Update the energy produced

        Parameters
        ----------
        fission_rates : numpy.ndarray
            fission reaction rate for each isotope in the specified
            material. Should be ordered corresponding to initial
            ``rate_index`` used in :meth:`prepare`
        mat_index : int
            Index for the specific material in the list of all burnable
            materials.
        """

    @property
    def nuclides(self):
        """List of nuclides with requested reaction rates"""
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        check_type("nuclides", nuclides, list, str)
        self._nuclides = nuclides
