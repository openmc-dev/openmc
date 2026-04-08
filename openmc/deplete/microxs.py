"""MicroXS module

A class for storing microscopic cross section data that can be used with the
IndependentOperator class for depletion.
"""

from __future__ import annotations
from collections.abc import Sequence
import shutil
from tempfile import TemporaryDirectory
from typing import Union, TypeAlias, Self

import h5py
import pandas as pd
import numpy as np

from openmc.checkvalue import check_type, check_value, check_iterable_type, PathLike
from openmc import StatePoint
from openmc.mgxs import GROUP_STRUCTURES
from openmc.data import REACTION_MT
import openmc
from .chain import Chain, REACTIONS, _get_chain
from .coupled_operator import _find_cross_sections, _get_nuclides_with_data
from ..utility_funcs import h5py_file_or_group
import openmc.lib
from openmc.mpi import comm

_valid_rxns = list(REACTIONS)
_valid_rxns.append('fission')
_valid_rxns.append('damage-energy')


# TODO: Replace with type statement when support is Python 3.12+
DomainTypes: TypeAlias = Union[
    Sequence[openmc.Material],
    Sequence[openmc.Cell],
    Sequence[openmc.Universe],
    openmc.MeshBase,
    openmc.Filter
]


def get_microxs_and_flux(
    model: openmc.Model,
    domains: DomainTypes,
    nuclides: Sequence[str] | None = None,
    reactions: Sequence[str] | None = None,
    energies: Sequence[float] | str | None = None,
    reaction_rate_mode: str = 'direct',
    chain_file: PathLike | Chain | None = None,
    path_statepoint: PathLike | None = None,
    path_input: PathLike | None = None,
    run_kwargs=None,
    reaction_rate_opts: dict | None = None,
) -> tuple[list[np.ndarray], list[MicroXS]]:
    """Generate microscopic cross sections and fluxes for multiple domains.

    This function runs a neutron transport solve to obtain the flux and reaction
    rates in the specified domains and computes multigroup microscopic cross
    sections that can be used in depletion calculations with the
    :class:`~openmc.deplete.IndependentOperator` class.

    .. versionadded:: 0.14.0

    .. versionchanged:: 0.15.3
        Added `reaction_rate_mode`, `path_statepoint`, `path_input` arguments.

    Parameters
    ----------
    model : openmc.Model
        OpenMC model object. Must contain geometry, materials, and settings.
    domains : list of openmc.Material or openmc.Cell or openmc.Universe, or openmc.MeshBase, or openmc.Filter
        Domains in which to tally reaction rates, or a spatial tally filter.
    nuclides : list of str
        Nuclides to get cross sections for. If not specified, all burnable
        nuclides from the depletion chain file are used.
    reactions : list of str
        Reactions to get cross sections for. If not specified, all neutron
        reactions listed in the depletion chain file are used.
    energies : iterable of float or str
        Energy group boundaries in [eV] or the name of the group structure.
        If left as None energies will default to [0.0, 100e6]
    reaction_rate_mode : {"direct", "flux"}, optional
        The "direct" method tallies reaction rates directly (per energy
        group). The "flux" method tallies a multigroup flux spectrum and then
        collapses reaction rates after a transport solve. When
        `reaction_rate_opts` is provided with `reaction_rate_mode='flux'`, the
        specified nuclide/reaction pairs are tallied directly and those values
        override the flux-collapsed values.
    chain_file : PathLike or Chain, optional
        Path to the depletion chain XML file or an instance of
        openmc.deplete.Chain. Used to determine cross sections for materials not
        present in the inital composition. Defaults to
        ``openmc.config['chain_file']``.
    path_statepoint : path-like, optional
        Path to write the statepoint file from the neutron transport solve to.
        By default, The statepoint file is written to a temporary directory and
        is not kept.
    path_input : path-like, optional
        Path to write the model XML file from the neutron transport solve to.
        By default, the model XML file is written to a temporary directory and
        not kept.
    run_kwargs : dict, optional
        Keyword arguments passed to :meth:`openmc.Model.run`
    reaction_rate_opts : dict, optional
        When `reaction_rate_mode="flux"`, allows selecting a subset of
        nuclide/reaction pairs to be computed via direct reaction-rate tallies
        (per energy group). Supported keys: "nuclides", "reactions".

    Returns
    -------
    list of numpy.ndarray
        Flux in each group in [n-cm/src] for each domain
    list of MicroXS
        Cross section data in [b] for each domain

    See Also
    --------
    openmc.deplete.IndependentOperator

    """
    check_value('reaction_rate_mode', reaction_rate_mode, {'direct', 'flux'})

    # Save any original tallies on the model
    original_tallies = list(model.tallies)

    # Determine what reactions and nuclides are available in chain
    chain = _get_chain(chain_file)
    if reactions is None:
        reactions = chain.reactions
    if not nuclides:
        cross_sections = _find_cross_sections(model)
        nuclides_with_data = _get_nuclides_with_data(cross_sections)
        nuclides = [nuc.name for nuc in chain.nuclides
                    if nuc.name in nuclides_with_data]

    # Set up the reaction rate and flux tallies
    if energies is None:
        energies = [0.0, 100.0e6]
    if isinstance(energies, str):
        energy_filter = openmc.EnergyFilter.from_group_structure(energies)
    else:
        energy_filter = openmc.EnergyFilter(energies)

    if isinstance(domains, openmc.Filter):
        domain_filter = domains
    elif isinstance(domains, openmc.MeshBase):
        domain_filter = openmc.MeshFilter(domains)
    elif isinstance(domains[0], openmc.Material):
        domain_filter = openmc.MaterialFilter(domains)
    elif isinstance(domains[0], openmc.Cell):
        domain_filter = openmc.CellFilter(domains)
    elif isinstance(domains[0], openmc.Universe):
        domain_filter = openmc.UniverseFilter(domains)
    else:
        raise ValueError(f"Unsupported domain type: {type(domains[0])}")

    flux_tally = openmc.Tally(name='MicroXS flux')
    flux_tally.filters = [domain_filter, energy_filter]
    flux_tally.scores = ['flux']
    model.tallies = [flux_tally]

    # Prepare reaction-rate tally for 'direct' or subset for 'flux' with opts
    rr_tally = None
    rr_nuclides: list[str] = []
    rr_reactions: list[str] = []
    if reaction_rate_mode == 'direct':
        rr_nuclides = list(nuclides)
        rr_reactions = list(reactions)
    elif reaction_rate_mode == 'flux' and reaction_rate_opts:
        opts = reaction_rate_opts or {}
        rr_nuclides = list(opts.get('nuclides', []))
        rr_reactions = list(opts.get('reactions', []))
        # Keep only requested pairs within overall sets
        if rr_nuclides:
            rr_nuclides = [n for n in rr_nuclides if n in set(nuclides)]
        if rr_reactions:
            rr_reactions = [r for r in rr_reactions if r in set(reactions)]

    # Only construct tally if both lists are non-empty
    if rr_nuclides and rr_reactions:
        rr_tally = openmc.Tally(name='MicroXS RR')
        # Use 1-group energy filter for RR in flux mode
        if reaction_rate_mode == 'flux':
            rr_energy_filter = openmc.EnergyFilter(
                [energy_filter.values[0], energy_filter.values[-1]])
        else:
            rr_energy_filter = energy_filter
        rr_tally.filters = [domain_filter, rr_energy_filter]
        rr_tally.nuclides = rr_nuclides
        rr_tally.multiply_density = False
        rr_tally.scores = rr_reactions
        model.tallies.append(rr_tally)

    if openmc.lib.is_initialized:
        openmc.lib.finalize()

        if comm.rank == 0:
            model.export_to_model_xml()
        comm.barrier()
        # Reinitialize with tallies
        openmc.lib.init(intracomm=comm)

    with TemporaryDirectory() as temp_dir:
        # Indicate to run in temporary directory unless being executed through
        # openmc.lib, in which case we don't need to specify the cwd
        run_kwargs = dict(run_kwargs) if run_kwargs else {}
        if not openmc.lib.is_initialized:
            run_kwargs.setdefault('cwd', temp_dir)

        # Run transport simulation and synchronize
        statepoint_path = model.run(**run_kwargs)
        comm.barrier()

        if comm.rank == 0:
            # Move the statepoint file if it is being saved to a specific path
            if path_statepoint is not None:
                shutil.move(statepoint_path, path_statepoint)
                statepoint_path = path_statepoint

            # Export the model to path_input if provided
            if path_input is not None:
                model.export_to_model_xml(path_input)

        # Broadcast updated statepoint path to all ranks
        statepoint_path = comm.bcast(statepoint_path)

        # Read in tally results (on all ranks)
        with StatePoint(statepoint_path) as sp:
            if rr_tally is not None:
                rr_tally = sp.tallies[rr_tally.id]
                rr_tally._read_results()
            flux_tally = sp.tallies[flux_tally.id]
            flux_tally._read_results()

    # Get flux values and make energy groups last dimension
    flux = flux_tally.get_reshaped_data()  # (domains, groups, 1, 1)
    flux = np.moveaxis(flux, 1, -1)  # (domains, 1, 1, groups)

    # Create list where each item corresponds to one domain
    fluxes = list(flux.squeeze((1, 2)))

    # If we built a reaction-rate tally, compute microscopic cross sections
    if rr_tally is not None:
        # Get reaction rates
        reaction_rates = rr_tally.get_reshaped_data()  # (domains, groups, nuclides, reactions)

        # Make energy groups last dimension
        reaction_rates = np.moveaxis(reaction_rates, 1, -1)  # (domains, nuclides, reactions, groups)

        # If RR is 1-group, sum flux over groups
        if reaction_rate_mode == "flux":
            flux = flux.sum(axis=-1, keepdims=True)  # (domains, 1, 1, 1)

        # Divide RR by flux to get microscopic cross sections. The indexing
        # ensures that only non-zero flux values are used, and broadcasting is
        # applied to align the shapes of reaction_rates and flux for division.
        xs = np.zeros_like(reaction_rates)  # (domains, nuclides, reactions, groups)
        d, _, _, g = np.nonzero(flux)
        xs[d, ..., g] = reaction_rates[d, ..., g] / flux[d, :, :, g]

        # Create lists where each item corresponds to one domain
        direct_micros = [MicroXS(xs_i, rr_nuclides, rr_reactions) for xs_i in xs]

    # If using flux mode, compute flux-collapsed microscopic XS
    if reaction_rate_mode == 'flux':
        flux_micros = [MicroXS.from_multigroup_flux(
            energies=energies,
            multigroup_flux=flux_i,
            chain_file=chain_file,
            nuclides=nuclides,
            reactions=reactions
        ) for flux_i in fluxes]

    # Decide which micros to use and merge if needed
    if reaction_rate_mode == 'flux' and rr_tally is not None:
        micros = [m1.merge(m2) for m1, m2 in zip(flux_micros, direct_micros)]
    elif rr_tally is not None:
        micros = direct_micros
    else:
        micros = flux_micros

    # Reset tallies
    model.tallies = original_tallies

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
    def __init__(self, data: np.ndarray, nuclides: list[str], reactions: list[str]):
        # Validate inputs
        if len(data.shape) != 3:
            raise ValueError('Data array must be 3D.')
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

    @classmethod
    def from_multigroup_flux(
        cls,
        energies: Sequence[float] | str,
        multigroup_flux: Sequence[float],
        chain_file: PathLike | None = None,
        temperature: float = 293.6,
        nuclides: Sequence[str] | None = None,
        reactions: Sequence[str] | None = None,
        **init_kwargs: dict,
    ) -> MicroXS:
        """Generated microscopic cross sections from a known flux.

        The size of the MicroXS matrix depends on the chain file and cross
        sections available. MicroXS entry will be 0 if the nuclide cross section
        is not found.

        It is recommended to make repeated calls to this method within a context
        manager using the :class:`openmc.lib.TemporarySession` class to avoid
        re-initializing OpenMC and loading cross sections each time.

        .. versionadded:: 0.15.0

        Parameters
        ----------
        energies : iterable of float or str
            Energy group boundaries in [eV] or the name of the group structure
        multigroup_flux : iterable of float
            Energy-dependent multigroup flux values
        chain_file : PathLike or Chain, optional
            Path to the depletion chain XML file or an instance of
            openmc.deplete.Chain. Defaults to ``openmc.config['chain_file']``.
        temperature : int, optional
            Temperature for cross section evaluation in [K].
        nuclides : list of str, optional
            Nuclides to get cross sections for. If not specified, all burnable
            nuclides from the depletion chain file are used.
        reactions : list of str, optional
            Reactions to get cross sections for. If not specified, all neutron
            reactions listed in the depletion chain file are used.
        **init_kwargs : dict
            Keyword arguments passed to :func:`openmc.lib.init`

        Returns
        -------
        MicroXS
        """

        check_type("temperature", temperature, (int, float))
        # if energy is string then use group structure of that name
        if isinstance(energies, str):
            energies = GROUP_STRUCTURES[energies]
        else:
            # if user inputs energies check they are ascending (low to high) as
            # some depletion codes use high energy to low energy.
            if not np.all(np.diff(energies) > 0):
                raise ValueError('Energy group boundaries must be in ascending order')

        # check dimension consistency
        if len(multigroup_flux) != len(energies) - 1:
            raise ValueError('Length of flux array should be len(energies)-1')

        chain = _get_chain(chain_file)
        cross_sections = _find_cross_sections(model=None)
        nuclides_with_data = _get_nuclides_with_data(cross_sections)

        # If no nuclides were specified, default to all nuclides from the chain
        if not nuclides:
            nuclides = chain.nuclides
            nuclides = [nuc.name for nuc in nuclides]

        # Get reaction MT values. If no reactions specified, default to the
        # reactions available in the chain file
        if reactions is None:
            reactions = chain.reactions
        mts = [REACTION_MT[name] for name in reactions]

        # Create 3D array for microscopic cross sections
        microxs_arr = np.zeros((len(nuclides), len(mts), 1))

        # If flux is zero, safely return zero cross sections
        multigroup_flux = np.array(multigroup_flux)
        if (flux_sum := multigroup_flux.sum()) == 0.0:
            return cls(microxs_arr, nuclides, reactions)

        # Normalize multigroup flux
        multigroup_flux /= flux_sum

        # Compute microscopic cross sections within a temporary session
        with openmc.lib.TemporarySession(**init_kwargs):
            # For each nuclide and reaction, compute the flux-averaged xs
            for nuc_index, nuc in enumerate(nuclides):
                if nuc not in nuclides_with_data:
                    continue
                lib_nuc = openmc.lib.load_nuclide(nuc)
                for mt_index, mt in enumerate(mts):
                    microxs_arr[nuc_index, mt_index, 0] = lib_nuc.collapse_rate(
                        mt, temperature, energies, multigroup_flux
                    )

        return cls(microxs_arr, nuclides, reactions)

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
        kwargs.setdefault('float_precision', 'round_trip')

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

    def to_hdf5(self, group_or_filename: h5py.Group | PathLike, **kwargs):
        """Export microscopic cross section data to HDF5 format

        Parameters
        ----------
        group_or_filename : h5py.Group or path-like
            HDF5 group or filename to write to
        kwargs : dict, optional
            Keyword arguments to pass to :meth:`h5py.Group.create_dataset`.
            Defaults to {'compression': 'lzf'}.

        """
        kwargs.setdefault('compression', 'lzf')

        with h5py_file_or_group(group_or_filename, 'w') as group:
            # Store cross section data as 3D dataset
            group.create_dataset('data', data=self.data, **kwargs)

            # Store metadata as datasets using string encoding
            group.create_dataset('nuclides', data=np.array(self.nuclides, dtype='S'))
            group.create_dataset('reactions', data=np.array(self.reactions, dtype='S'))

    @classmethod
    def from_hdf5(cls, group_or_filename: h5py.Group | PathLike) -> Self:
        """Load data from an HDF5 file

        Parameters
        ----------
        group_or_filename : h5py.Group or str or PathLike
            HDF5 group or path to HDF5 file. If given as an h5py.Group, the
            data is read from that group. If given as a string, it is assumed
            to be the filename for the HDF5 file.

        Returns
        -------
        MicroXS
        """

        with h5py_file_or_group(group_or_filename, 'r') as group:
            # Read data from HDF5 group
            data = group['data'][:]
            nuclides = [nuc.decode('utf-8') for nuc in group['nuclides'][:]]
            reactions = [rxn.decode('utf-8') for rxn in group['reactions'][:]]

        return cls(data, nuclides, reactions)

    def merge(self, other: Self, prefer: str = 'other') -> Self:
        """Merge two MicroXS objects by taking the union of nuclides/reactions.

        If the two objects contain overlapping nuclide/reaction entries, values
        from `other` will overwrite values from `self` when `prefer='other'`.
        When `prefer='self'`, values from `self` are retained for overlapping
        entries, and values from `other` are used only for non-overlapping
        entries.

        Parameters
        ----------
        other : MicroXS
            Other MicroXS instance to merge with this one.
        prefer : {"other", "self"}
            Which instance's data should take precedence on overlap.

        Returns
        -------
        MicroXS
            New instance containing the merged data.
        """
        check_value('prefer', prefer, {'other', 'self'})

        # Require same number of energy groups
        if self.data.shape[2] != other.data.shape[2]:
            raise ValueError(
                'Cannot merge MicroXS with different number of energy groups: '
                f"{self.data.shape[2]} vs {other.data.shape[2]}. Ensure that "
                'both were generated with consistent group structures and '
                'treatments (e.g., both multigroup or both collapsed).'
            )

        # Build unified axes preserving order (self first, then other's new)
        new_nuclides = list(self.nuclides)
        for nuc in other.nuclides:
            if nuc not in self._index_nuc:
                new_nuclides.append(nuc)
        new_reactions = list(self.reactions)
        for rx in other.reactions:
            if rx not in self._index_rx:
                new_reactions.append(rx)

        # Allocate and fill from self (self's nuclides/reactions map to the
        # first indices of new_nuclides/new_reactions by construction)
        groups = self.data.shape[2]
        data = np.zeros((len(new_nuclides), len(new_reactions), groups))
        idx_n = {nuc: i for i, nuc in enumerate(new_nuclides)}
        idx_r = {rx: i for i, rx in enumerate(new_reactions)}

        n_self = len(self.nuclides)
        r_self = len(self.reactions)
        data[:n_self, :r_self] = self.data

        # Build destination index arrays for other's nuclides/reactions
        dst_n = np.array([idx_n[nuc] for nuc in other.nuclides])
        dst_r = np.array([idx_r[rx] for rx in other.reactions])

        # Copy from other, respecting precedence
        if prefer == 'other':
            data[np.ix_(dst_n, dst_r)] = other.data
        else:
            # Copy only entries where nuc or rx is absent from self
            nuc_is_new = np.array(
                [nuc not in self._index_nuc for nuc in other.nuclides])
            rx_is_new = np.array(
                [rx not in self._index_rx for rx in other.reactions])
            mask = nuc_is_new[:, np.newaxis] | rx_is_new[np.newaxis, :]
            src_i, src_j = np.where(mask)
            if src_i.size:
                data[dst_n[src_i], dst_r[src_j]] = other.data[src_i, src_j]

        return MicroXS(data, new_nuclides, new_reactions)


def write_microxs_hdf5(
    micros: Sequence[MicroXS],
    filename: PathLike,
    names: Sequence[str] | None = None,
    **kwargs
):
    """Write multiple MicroXS objects to an HDF5 file

    Parameters
    ----------
    micros : list of MicroXS
        List of MicroXS objects
    filename : PathLike
        Output HDF5 filename
    names : list of str, optional
        Names for each MicroXS object. If None, uses 'domain_0', 'domain_1',
        etc.
    **kwargs
        Additional keyword arguments passed to :meth:`h5py.Group.create_dataset`
    """
    if names is None:
        names = [f'domain_{i}' for i in range(len(micros))]

    # Open file once and write all domains using group interface
    with h5py.File(filename, 'w') as f:
        for microxs, name in zip(micros, names):
            group = f.create_group(name)
            microxs.to_hdf5(group, **kwargs)


def read_microxs_hdf5(filename: PathLike) -> dict[str, MicroXS]:
    """Read multiple MicroXS objects from an HDF5 file

    Parameters
    ----------
    filename : path-like
        HDF5 filename

    Returns
    -------
    dict
        Dictionary mapping domain names to MicroXS objects
    """
    with h5py.File(filename, 'r') as f:
        return {name: MicroXS.from_hdf5(group) for name, group in f.items()}
