"""MicroXS module

A class for storing microscopic cross section data that can be used with the
IndependentOperator class for depletion.
"""

import tempfile
from typing import List, Iterable, Optional

import pandas as pd
import numpy as np

from openmc.checkvalue import check_type, check_value, check_iterable_type, PathLike
from openmc.exceptions import DataError
from openmc import StatePoint
import openmc
from .chain import Chain, REACTIONS
from .coupled_operator import _find_cross_sections, _get_nuclides_with_data

_valid_rxns = list(REACTIONS)
_valid_rxns.append('fission')


def get_microxs_and_flux(
        model: openmc.Model,
        domain,
        nuclides: Optional[Iterable[str]] = None,
        reactions: Optional[Iterable[str]] = None,
        chain_file: Optional[PathLike] =None,
        energy_bounds: Optional[Iterable[float]] = (0, 20e6),
        run_kwargs=None
    ):
    """Generate a microscopic cross sections and flux from a Model

    Parameters
    ----------
    model : openmc.Model
        OpenMC model object. Must contain geometry, materials, and settings.
    domain : openmc.Material or openmc.Cell or openmc.Universe
        Domain in which to tally reaction rates.
    nuclides : list of str
        Nuclides to get cross sections for. If not specified, all burnable
        nuclides from the depletion chain file are used.
    reactions : list of str
        Reactions to get cross sections for. If not specified, all neutron
        reactions listed in the depletion chain file are used.
    chain_file : str, optional
        Path to the depletion chain XML file that will be used in depletion
        simulation. Used to determine cross sections for materials not
        present in the inital composition. Defaults to
        ``openmc.config['chain_file']``.
    energy_bound : 2-tuple of float, optional
        Bounds for the energy group.
    run_kwargs : dict, optional
        Keyword arguments passed to :meth:`openmc.Model.run`

    Returns
    -------
    float
        One-group flux in [n-cm/src]
    MicroXS
        Cross section data in [b]

    """
    # Save any original tallies on the model
    original_tallies = model.tallies

    # Determine what reactions and nuclides are available in chain
    if chain_file is None:
        chain_file = openmc.config.get('chain_file')
        if chain_file is None:
            raise DataError(
                "No depletion chain specified and could not find depletion "
                "chain in openmc.config['chain_file']"
            )
    chain = Chain.from_xml(chain_file)
    if reactions is None:
        reactions = chain.reactions
    if not nuclides:
        cross_sections = _find_cross_sections(model)
        nuclides_with_data = _get_nuclides_with_data(cross_sections)
        nuclides = [nuc.name for nuc in chain.nuclides
                    if nuc.name in nuclides_with_data]

    # Set up the reaction rate and flux tallies
    energy_filter = openmc.EnergyFilter(energy_bounds)
    if isinstance(domain, openmc.Material):
        domain_filter = openmc.MaterialFilter([domain])
    elif isinstance(domain, openmc.Cell):
        domain_filter = openmc.CellFilter([domain])
    elif isinstance(domain, openmc.Universe):
        domain_filter = openmc.UniverseFilter([domain])
    else:
        raise ValueError(f"Unsupported domain type: {type(domain)}")

    rr_tally = openmc.Tally(name='MicroXS RR')
    rr_tally.filters = [domain_filter, energy_filter]
    rr_tally.nuclides = nuclides
    rr_tally.multiply_density = False
    rr_tally.scores = reactions

    flux_tally = openmc.Tally(name='MicroXS flux')
    flux_tally.filters = [domain_filter, energy_filter]
    flux_tally.scores = ['flux']
    model.tallies = openmc.Tallies([rr_tally, flux_tally])

    # create temporary run
    with tempfile.TemporaryDirectory() as temp_dir:
        if run_kwargs is None:
            run_kwargs = {}
        run_kwargs.setdefault('cwd', temp_dir)
        statepoint_path = model.run(**run_kwargs)

        with StatePoint(statepoint_path) as sp:
            rr_tally = sp.tallies[rr_tally.id]
            rr_tally._read_results()
            flux_tally = sp.tallies[flux_tally.id]
            flux_tally._read_results()

    # Get reaction rates and flux values
    reaction_rates = rr_tally.mean.sum(axis=0)  # (nuclides, reactions)
    flux = flux_tally.mean[0, 0, 0]

    # Divide RR by flux to get microscopic cross sections
    xs = reaction_rates / flux  # (nuclides, reactions)

    # Reset tallies
    model.tallies = original_tallies

    return flux, MicroXS(xs, nuclides, reactions)


class MicroXS:
    """Microscopic cross section data for use in transport-independent depletion.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    data : numpy.ndarray of floats
        Array containing one-group microscopic cross section values for
        each nuclide and reaction. Cross section values are assumed to be
        in [b].
    nuclides : list of str
        List of nuclide symbols for that have data for at least one
        reaction.
    reactions : list of str
        List of reactions. All reactions must match those in
        :data:`openmc.deplete.chain.REACTIONS`

    """
    def __init__(self, data: np.ndarray, nuclides: List[str], reactions: List[str]):
        # Validate inputs
        if data.shape != (len(nuclides), len(reactions)):
            raise ValueError(
                f'Nuclides list of length {len(nuclides)} and '
                f'reactions array of length {len(reactions)} do not '
                f'match dimensions of data array of shape {data.shape}')
        check_iterable_type('nuclides', nuclides, str)
        check_iterable_type('reactions', reactions, str)
        check_type('data', data, np.ndarray, expected_iter_type=float)
        for reaction in reactions:
            check_value('reactions', reaction, _valid_rxns)

        self.data = data
        self.nuclides = nuclides
        self.reactions = reactions

    @classmethod
    def from_csv(cls, csv_file, **kwargs):
        """Load data from a comma-separated values (csv) file.

        Parameters
        ----------
        csv_file : str
            Relative path to csv-file containing microscopic cross section
            data. Cross section values are assumed to be in [b]
        **kwargs : dict
            Keyword arguments to pass to :func:`pandas.read_csv()`.

        Returns
        -------
        MicroXS

        """
        if 'float_precision' not in kwargs:
            kwargs['float_precision'] = 'round_trip'

        df = pd.read_csv(csv_file, index_col=0, **kwargs)
        data = df.values
        nuclides = list(df.index)
        reactions = list(df)
        return cls(data, nuclides, reactions)

    def __getitem__(self, index):
        nuc, rx = index
        i_nuc = self.nuclides.index(nuc)
        i_rx = self.reactions.index(rx)
        return self.data[i_nuc, i_rx]

    def to_csv(self, *args, **kwargs):
        """Write data to a comma-separated values (csv) file

        Parameters
        ----------
        *args
            Positional arguments passed to :meth:`pandas.DataFrame.to_csv`
        **kwargs
            Keyword arguments passed to :meth:`pandas.DataFrame.to_csv`

        """
        index = pd.Index(self.nuclides, name='nuclide')
        df = pd.DataFrame(self.data, columns=self.reactions, index=index)
        df.to_csv(*args, **kwargs)

