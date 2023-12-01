"""MicroXS module

A class for storing microscopic cross section data that can be used with the
IndependentOperator class for depletion.
"""

from __future__ import annotations
import tempfile
from typing import List, Tuple, Iterable, Optional, Union

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
        domains,
        nuclides: Optional[Iterable[str]] = None,
        reactions: Optional[Iterable[str]] = None,
        energies: Optional[Union[Iterable[float], str]] = None,
        chain_file: Optional[PathLike] = None,
        run_kwargs=None
    ) -> Tuple[List[np.ndarray], List[MicroXS]]:
    """Generate a microscopic cross sections and flux from a Model

    .. versionadded:: 0.14.0

    Parameters
    ----------
    model : openmc.Model
        OpenMC model object. Must contain geometry, materials, and settings.
    domains : list of openmc.Material or openmc.Cell or openmc.Universe, or openmc.MeshBase
        Domains in which to tally reaction rates.
    nuclides : list of str
        Nuclides to get cross sections for. If not specified, all burnable
        nuclides from the depletion chain file are used.
    reactions : list of str
        Reactions to get cross sections for. If not specified, all neutron
        reactions listed in the depletion chain file are used.
    energies : iterable of float or str
        Energy group boundaries in [eV] or the name of the group structure
    chain_file : str, optional
        Path to the depletion chain XML file that will be used in depletion
        simulation. Used to determine cross sections for materials not
        present in the inital composition. Defaults to
        ``openmc.config['chain_file']``.
    run_kwargs : dict, optional
        Keyword arguments passed to :meth:`openmc.Model.run`

    Returns
    -------
    list of numpy.ndarray
        Flux in each group in [n-cm/src] for each domain
    list of MicroXS
        Cross section data in [b] for each domain

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
    if energies is None:
        energy_filter = openmc.EnergyFilter([0.0, 100.0e6])
    elif isinstance(energies, str):
        energy_filter = openmc.EnergyFilter.from_group_structure(energies)
    else:
        energy_filter = openmc.EnergyFilter(energies)
    if isinstance(domains, openmc.MeshBase):
        domain_filter = openmc.MeshFilter(domains)
    elif isinstance(domains[0], openmc.Material):
        domain_filter = openmc.MaterialFilter(domains)
    elif isinstance(domains[0], openmc.Cell):
        domain_filter = openmc.CellFilter(domains)
    elif isinstance(domains[0], openmc.Universe):
        domain_filter = openmc.UniverseFilter(domains)
    else:
        raise ValueError(f"Unsupported domain type: {type(domains[0])}")

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
        else:
            run_kwargs = dict(run_kwargs)
        run_kwargs.setdefault('cwd', temp_dir)
        statepoint_path = model.run(**run_kwargs)

        with StatePoint(statepoint_path) as sp:
            rr_tally = sp.tallies[rr_tally.id]
            rr_tally._read_results()
            flux_tally = sp.tallies[flux_tally.id]
            flux_tally._read_results()

    # Get reaction rates and flux values
    reaction_rates = rr_tally.get_reshaped_data()  # (domains, groups, nuclides, reactions)
    flux = flux_tally.get_reshaped_data()  # (domains, groups, 1, 1)

    # Make energy groups last dimension
    reaction_rates = np.moveaxis(reaction_rates, 1, -1)  # (domains, nuclides, reactions, groups)
    flux = np.moveaxis(flux, 1, -1)  # (domains, 1, 1, groups)

    # Divide RR by flux to get microscopic cross sections
    xs = np.empty_like(reaction_rates) # (domains, nuclides, reactions, groups)
    d, _, _, g = np.nonzero(flux)
    xs[d, ..., g] = reaction_rates[d, ..., g] / flux[d, :, :, g]

    # Reset tallies
    model.tallies = original_tallies

    # Create lists where each item corresponds to one domain
    fluxes = list(flux.squeeze((1, 2)))
    micros = [MicroXS(xs_i, nuclides, reactions) for xs_i in xs]
    return fluxes, micros


class MicroXS:
    """Microscopic cross section data for use in transport-independent depletion.

    .. versionadded:: 0.13.1

    .. versionchanged:: 0.14.0
        Class was heavily refactored and no longer subclasses pandas.DataFrame.

    Parameters
    ----------
    data : numpy.ndarray of floats
        3D array containing microscopic cross section values for each
        nuclide, reaction, and energy group. Cross section values are assumed to
        be in [b], and indexed by [nuclide, reaction, energy group]
    nuclides : list of str
        List of nuclide symbols for that have data for at least one
        reaction.
    reactions : list of str
        List of reactions. All reactions must match those in
        :data:`openmc.deplete.chain.REACTIONS`

    """
    def __init__(self, data: np.ndarray, nuclides: List[str], reactions: List[str]):
        # Validate inputs
        if data.shape[:2] != (len(nuclides), len(reactions)):
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
        self._index_nuc = {nuc: i for i, nuc in enumerate(nuclides)}
        self._index_rx = {rx: i for i, rx in enumerate(reactions)}

    # TODO: Add a classmethod for generating MicroXS directly from cross section
    # data using a known flux spectrum

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

        df = pd.read_csv(csv_file, **kwargs)
        df.set_index(['nuclides', 'reactions', 'groups'], inplace=True)
        nuclides = list(df.index.unique(level='nuclides'))
        reactions = list(df.index.unique(level='reactions'))
        groups = list(df.index.unique(level='groups'))
        shape = (len(nuclides), len(reactions), len(groups))
        data = df.values.reshape(shape)
        return cls(data, nuclides, reactions)

    def __getitem__(self, index):
        nuc, rx = index
        i_nuc = self._index_nuc[nuc]
        i_rx = self._index_rx[rx]
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
        groups = self.data.shape[2]
        multi_index = pd.MultiIndex.from_product(
            [self.nuclides, self.reactions, range(1, groups + 1)],
            names=['nuclides', 'reactions', 'groups']
        )
        df = pd.DataFrame({'xs': self.data.flatten()}, index=multi_index)
        df.to_csv(*args, **kwargs)

