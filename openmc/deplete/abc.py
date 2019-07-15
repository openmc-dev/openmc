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

from numpy import zeros, nonzero

from openmc.data import DataLibrary
from openmc.checkvalue import check_type
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

    Attributes
    ----------
    dilute_initial : float
        Initial atom density to add for nuclides that are zero in initial
        condition to ensure they exist in the decay chain.  Only done for
        nuclides with reaction rates. Defaults to 1.0e3.

    """
    def __init__(self, chain_file=None, fission_q=None):
        self.dilute_initial = 1.0e3
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

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates. Ordered to be
        consistent with :class:`openmc.deplete.Operator`
    """

    def __init__(self):
        self._nuclides = None
        self._rate_tally = None
        self._results_cache = None

    @abstractmethod
    def generate_tallies(self, materials, scores):
        """Use the capi to build tallies needed for reaction rates"""

    @property
    def nuclides(self):
        """List of nuclides with requested reaction rates"""
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        check_type("nuclides", nuclides, list, str)
        self._nuclides = nuclides
        self._rate_tally.nuclides = nuclides

    def _reset_results_cache(self, nnucs, nreact):
        """Cache for results for a given material"""
        if (self._results_cache is None
                or self._results_cache.shape != (nnucs, nreact)):
            self._results_cache = zeros((nnucs, nreact))
        else:
            self._results_cache.fill(0.0)
        return self._results_cache

    @abstractmethod
    def get_material_rates(self, mat_id, nuc_index, react_index):
        """Return 2D array of [nuclide, reaction] reaction rates

        ``nuc_index`` and ``react_index`` are orderings of nuclides
        and reactions such that the ordering is consistent between
        reaction tallies and energy deposition tallies"""
        pass

    def divide_by_adens(self, number):
        """Normalize reaction rates by number of nuclides

        Acts on the current material examined by
        :meth:`get_material_rates`

        Parameters
        ----------
        energy : float
            Energy produced in this region [W]
        number : iterable of float
            Number density [#/b/cm] of each nuclide tracked in the calculation.
            Ordered identically to :attr:`nuclides`

        Returns
        -------
        results : :class:`numpy.ndarray`
            2D array ``[n_nuclides, n_rxns]`` of reaction rates normalized by
            the number of nuclides
        """

        mask = nonzero(number)
        results = self._results_cache
        for col in range(results.shape[1]):
            results[mask, col] /= number[mask]
        return results


class FissionEnergyHelper(ABC):
    """Abstract class for normalizing fission reactions to a given level

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates. Ordered to be
        consistent with :class:`openmc.deplete.Operator`
    """

    def __init__(self):
        self._nuclides = None
        self._fission_E = None

    @abstractmethod
    def prepare(self, chain_nucs, rate_index, materials):
        """Perform work needed to obtain fission energy per material

        ``chain_nucs`` is all nuclides tracked in the depletion chain,
        while ``rate_index`` should be a mapping from nuclide name
        to index in the reaction rate vector used in
        :meth:`get_fission_energy`.
        ``materials`` should be a list of all materials tracked
        on the operator to which this object is attached"""

    @abstractmethod
    def get_fission_energy(self, fission_rates, mat_index):
        """Return fission energy in this material given fission rates"""

    @property
    def nuclides(self):
        """List of nuclides with requested reaction rates"""
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        check_type("nuclides", nuclides, list, str)
        self._nuclides = nuclides
