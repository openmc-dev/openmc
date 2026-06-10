"""MicroXS module

A class for storing microscopic cross section data that can be used with the
IndependentOperator class for depletion.
"""

from __future__ import annotations
from collections.abc import Sequence
from contextlib import nullcontext
from dataclasses import dataclass
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
    openmc.Filter,
    Sequence[openmc.Filter]
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
    domains : list of openmc.Material or openmc.Cell or openmc.Universe, or openmc.MeshBase, or openmc.Filter, or list of openmc.Filter
        Domains in which to tally reaction rates, or a spatial tally filter.
        A list of filters can be provided to create one set of tallies per
        filter (e.g., one :class:`~openmc.MeshMaterialFilter` per mesh) that
        are all evaluated in a single transport solve. Results are
        concatenated across all filters in order.
    nuclides : list of str
        Nuclides to get cross sections for. If not specified, all burnable
        nuclides from the depletion chain file are used.
    reactions : list of str
        Reactions to get cross sections for. If not specified, all neutron
        reactions listed in the depletion chain file are used.
    energies : iterable of float or str
        Energy group boundaries in [eV] or the name of the group structure.
        If left as None, no energy filter is applied to the flux tally. When
        `reaction_rate_mode` is "direct", these boundaries define the output
        flux and microscopic cross section energy group structure. When
        `reaction_rate_mode` is "flux", these boundaries define the multigroup
        flux tally used to collapse continuous-energy cross sections; returned
        fluxes and microscopic cross sections are one-group.
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
        over one energy bin spanning the full `energies` range. Supported keys:
        "nuclides", "reactions". If "reactions" are specified without
        "nuclides", all selected nuclides are used.

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

    # Set up the reaction rate and flux tallies. When energies are omitted, no
    # energy filter is needed for the transport calculation. A one-group energy
    # range is still needed later if flux collapse is requested.
    collapse_energies = energies
    if energies is None:
        energy_filter = None
        collapse_energies = [0.0, 100.0e6]
    elif isinstance(energies, str):
        energy_filter = openmc.EnergyFilter.from_group_structure(energies)
    else:
        energy_filter = openmc.EnergyFilter(energies)

    # Build list of domain filters
    if isinstance(domains, openmc.Filter):
        domain_filters = [domains]
    elif isinstance(domains, openmc.MeshBase):
        domain_filters = [openmc.MeshFilter(domains)]
    elif isinstance(domains, Sequence) and len(domains) > 0 and \
            isinstance(domains[0], openmc.Filter):
        domain_filters = list(domains)
    elif isinstance(domains[0], openmc.Material):
        domain_filters = [openmc.MaterialFilter(domains)]
    elif isinstance(domains[0], openmc.Cell):
        domain_filters = [openmc.CellFilter(domains)]
    elif isinstance(domains[0], openmc.Universe):
        domain_filters = [openmc.UniverseFilter(domains)]
    else:
        raise ValueError(f"Unsupported domain type: {type(domains[0])}")

    # Prepare reaction-rate nuclides/reactions
    rr_nuclides: list[str] = []
    rr_reactions: list[str] = []
    if reaction_rate_mode == 'direct':
        rr_nuclides = list(nuclides)
        rr_reactions = list(reactions)
    elif reaction_rate_mode == 'flux' and reaction_rate_opts:
        opts = reaction_rate_opts or {}
        rr_reactions = list(opts.get('reactions', []))
        if rr_reactions:
            rr_nuclides = list(opts.get('nuclides', nuclides))
        else:
            rr_nuclides = list(opts.get('nuclides', []))
        # Keep only requested pairs within overall sets
        if rr_nuclides:
            rr_nuclides = [n for n in rr_nuclides if n in set(nuclides)]
        if rr_reactions:
            rr_reactions = [r for r in rr_reactions if r in set(reactions)]

    # Use 1-group energy filter for RR in flux mode
    has_rr = bool(rr_nuclides and rr_reactions)
    if has_rr and reaction_rate_mode == 'flux' and energy_filter is not None:
        rr_energy_filter = openmc.EnergyFilter(
            [energy_filter.values[0], energy_filter.values[-1]])
    else:
        rr_energy_filter = energy_filter

    # Create one flux tally (and optionally one RR tally) per domain filter.
    flux_tallies = []
    rr_tallies = []
    model.tallies = []
    for i, domain_filter in enumerate(domain_filters):
        flux_tally = openmc.Tally(name=f'MicroXS flux {i}')
        flux_tally.filters = [domain_filter]
        if energy_filter is not None:
            flux_tally.filters.append(energy_filter)
        flux_tally.scores = ['flux']
        model.tallies.append(flux_tally)
        flux_tallies.append(flux_tally)

        if has_rr:
            rr_tally = openmc.Tally(name=f'MicroXS RR {i}')
            rr_tally.filters = [domain_filter]
            if rr_energy_filter is not None:
                rr_tally.filters.append(rr_energy_filter)
            rr_tally.nuclides = rr_nuclides
            rr_tally.multiply_density = False
            rr_tally.scores = rr_reactions
            model.tallies.append(rr_tally)
            rr_tallies.append(rr_tally)

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
            for i in range(len(flux_tallies)):
                flux_tallies[i] = sp.tallies[flux_tallies[i].id]
                flux_tallies[i]._read_results()
                if rr_tallies:
                    rr_tallies[i] = sp.tallies[rr_tallies[i].id]
                    rr_tallies[i]._read_results()

    # Concatenate results across all domain filters
    fluxes = []
    all_flux_arrays = []
    for flux_tally in flux_tallies:
        # Get flux values and make energy groups last dimension
        flux = flux_tally.get_reshaped_data()
        if energy_filter is None:
            flux = flux[..., np.newaxis]  # (domains, 1, 1, groups)
        else:
            # (domains, groups, 1, 1) -> (domains, 1, 1, groups)
            flux = np.moveaxis(flux, 1, -1)
        all_flux_arrays.append(flux)
        fluxes.extend(flux.squeeze((1, 2)))

    # If we built reaction-rate tallies, compute microscopic cross sections
    if rr_tallies:
        direct_micros = []
        for flux_arr, rr_tally in zip(all_flux_arrays, rr_tallies):
            flux = flux_arr
            # Get reaction rates and make energy groups last dimension
            reaction_rates = rr_tally.get_reshaped_data()
            if rr_energy_filter is None:
                # (domains, nuclides, reactions) ->
                # (domains, nuclides, reactions, groups)
                reaction_rates = reaction_rates[..., np.newaxis]
            else:
                # (domains, groups, nuclides, reactions) ->
                # (domains, nuclides, reactions, groups)
                reaction_rates = np.moveaxis(reaction_rates, 1, -1)

            # If RR is 1-group, sum flux over groups
            if reaction_rate_mode == "flux":
                flux = flux.sum(axis=-1, keepdims=True)

            xs = np.zeros_like(reaction_rates)
            d, _, _, g = np.nonzero(flux)
            xs[d, ..., g] = reaction_rates[d, ..., g] / flux[d, :, :, g]
            direct_micros.extend(
                MicroXS(xs_i, rr_nuclides, rr_reactions) for xs_i in xs)

    if reaction_rate_mode == 'flux':
        # Resolve the library from the model (from_multigroup_flux defaults to config)
        cross_sections = _find_cross_sections(model)
        # Collapse all domains against one table, built once
        flux_micros = MicroXS.from_multigroup_flux(
            collapse_energies, fluxes, chain_file=chain, nuclides=nuclides,
            reactions=reactions, cross_sections=cross_sections)

        # We need to return one-group fluxes to match the microscopic cross
        # sections, which are always one-group by virtue of the collapse
        fluxes = [flux.sum(keepdims=True) for flux in fluxes]

    # Decide which micros to use and merge if needed
    if reaction_rate_mode == 'flux' and rr_tallies:
        micros = [m1.merge(m2) for m1, m2 in zip(flux_micros, direct_micros)]
    elif rr_tallies:
        micros = direct_micros
    else:
        micros = flux_micros

    # Reset tallies
    model.tallies = original_tallies

    return fluxes, micros


@dataclass
class _SparseXSTable:
    """Sparse storage of group cross sections for vectorized flux collapse.

    Stores only the non-zero ``(nuclide, reaction)`` pairs. ``xs_matrix`` has
    shape ``(nnz, n_groups)``; ``nuc_indices`` and ``rxn_indices`` map each row
    to its position in the dense ``(n_nuclides, n_reactions)`` result.

    Parameters
    ----------
    nuclides : list of str
        Nuclide names defining the result's nuclide axis.
    reactions : list of str
        Reaction names defining the result's reaction axis.
    n_groups : int
        Number of energy groups.
    xs_matrix : numpy.ndarray
        Group cross sections for the stored rows, shape ``(nnz, n_groups)``.
    nuc_indices : numpy.ndarray
        Row-to-nuclide index map, shape ``(nnz,)``.
    rxn_indices : numpy.ndarray
        Row-to-reaction index map, shape ``(nnz,)``.

    Notes
    -----
    Provider-agnostic: rows may be filled from continuous-energy data
    (:func:`_build_xs_table_ce`) or any other group cross section source.

    """
    nuclides: list[str]
    reactions: list[str]
    n_groups: int
    xs_matrix: np.ndarray
    nuc_indices: np.ndarray
    rxn_indices: np.ndarray

    def collapse(self, phi_norm: np.ndarray) -> np.ndarray:
        """Collapse the table to one-group cross sections.

        Parameters
        ----------
        phi_norm : numpy.ndarray
            Normalized group flux, shape ``(n_groups,)``. Should sum to 1.

        Returns
        -------
        numpy.ndarray
            Dense ``(n_nuclides, n_reactions)`` one-group cross sections.

        """
        if len(phi_norm) != self.n_groups:
            raise ValueError(
                f'phi_norm has length {len(phi_norm)} but table expects '
                f'{self.n_groups} groups')
        collapsed_sparse = self.xs_matrix @ phi_norm
        result = np.zeros((len(self.nuclides), len(self.reactions)))
        result[self.nuc_indices, self.rxn_indices] = collapsed_sparse
        return result


def _build_xs_table_ce(
    nuclides: Sequence[str],
    reactions: Sequence[str],
    energies: Sequence[float],
    temperature: float,
    nuclides_with_data: set,
    cross_sections=None,
    **init_kwargs,
) -> _SparseXSTable:
    """Build a sparse group cross section table from continuous-energy data.

    For each requested ``(nuclide, reaction)`` the group-averaged cross sections
    are computed once via the :meth:`openmc.lib.Nuclide.group_xs` C API and
    stored as one sparse row, so every domain can later be collapsed with a
    single matrix-vector product. Rows that are entirely zero (reaction absent,
    or threshold above the whole group structure) are skipped. Cross sections
    are evaluated inside a single :class:`openmc.lib.TemporarySession`, loading
    each nuclide once.

    Parameters
    ----------
    nuclides : sequence of str
        Nuclide names defining the result's nuclide axis.
    reactions : sequence of str
        Reaction names defining the result's reaction axis.
    energies : sequence of float
        Energy group boundaries in [eV], ascending, length ``n_groups + 1``.
    temperature : float
        Temperature in [K] for cross section evaluation.
    nuclides_with_data : set
        Nuclides available in the cross section library; others are skipped.
    cross_sections : PathLike, optional
        Cross section library for the session, so it matches the one
        ``nuclides_with_data`` was resolved from. Defaults to ``openmc.config``.
    **init_kwargs : dict
        Keyword arguments passed to :func:`openmc.lib.init` via the temporary
        session.

    Returns
    -------
    _SparseXSTable

    """
    mts = [REACTION_MT[name] for name in reactions]
    energies = np.asarray(energies, dtype=float)
    n_groups = len(energies) - 1

    rows = []
    nuc_idx_list = []
    rxn_idx_list = []
    # Load against the same library nuclides_with_data was resolved from
    library = (openmc.config.patch('cross_sections', cross_sections)
               if cross_sections is not None else nullcontext())
    with library, openmc.lib.TemporarySession(**init_kwargs):
        for nuc_idx, nuc in enumerate(nuclides):
            if nuc not in nuclides_with_data:
                continue
            lib_nuc = openmc.lib.load_nuclide(nuc)
            # Index by reaction, not MT, so fission/(n,fission) stay separate
            for rxn_idx, mt in enumerate(mts):
                xs_g = lib_nuc.group_xs(mt, temperature, energies)
                if xs_g.any():
                    rows.append(xs_g)
                    nuc_idx_list.append(nuc_idx)
                    rxn_idx_list.append(rxn_idx)

    if rows:
        xs_matrix = np.vstack(rows)
    else:
        xs_matrix = np.empty((0, n_groups))

    return _SparseXSTable(
        nuclides=list(nuclides),
        reactions=list(reactions),
        n_groups=n_groups,
        xs_matrix=xs_matrix,
        nuc_indices=np.array(nuc_idx_list, dtype=np.int32),
        rxn_indices=np.array(rxn_idx_list, dtype=np.int32),
    )


def _collapse_fluxes(
    table: _SparseXSTable,
    fluxes: Sequence[np.ndarray],
    nuclides: Sequence[str],
    reactions: Sequence[str],
) -> list[MicroXS]:
    """Collapse each domain's multigroup flux against a built XS table.

    Each flux is validated (finite, non-negative) and normalized to sum 1 before
    being collapsed; a zero-sum flux yields an all-zero :class:`MicroXS`.

    Parameters
    ----------
    table : _SparseXSTable
        Pre-built group cross section table.
    fluxes : sequence of numpy.ndarray
        One ``(n_groups,)`` group-flux vector per domain.
    nuclides : sequence of str
        Nuclide axis of the output MicroXS objects.
    reactions : sequence of str
        Reaction axis of the output MicroXS objects.

    Returns
    -------
    list of MicroXS
        One per domain, each of shape ``(n_nuclides, n_reactions, 1)``.

    """
    micros = []
    for flux in fluxes:
        flux = np.asarray(flux, dtype=float)
        if not np.isfinite(flux).all():
            raise ValueError('Multigroup flux contains non-finite values')
        if (flux < 0).any():
            raise ValueError('Multigroup flux contains negative values')
        flux_sum = flux.sum()
        if flux_sum == 0.0:
            micros.append(MicroXS(
                np.zeros((len(nuclides), len(reactions), 1)),
                nuclides, reactions))
            continue
        collapsed = table.collapse(flux / flux_sum)
        micros.append(MicroXS(collapsed[:, :, np.newaxis],
                              nuclides, reactions))
    return micros


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
        multigroup_flux: Sequence[float] | Sequence[Sequence[float]],
        chain_file: PathLike | None = None,
        temperature: float = 293.6,
        nuclides: Sequence[str] | None = None,
        reactions: Sequence[str] | None = None,
        *,
        cross_sections: PathLike | None = None,
        **init_kwargs: dict,
    ) -> MicroXS | list[MicroXS]:
        """Generated microscopic cross sections from a known flux.

        The size of the MicroXS matrix depends on the chain file and cross
        sections available. MicroXS entry will be 0 if the nuclide cross section
        is not found.

        Multiple fluxes can be collapsed at once by passing a 2-D array (or a
        list of 1-D arrays); the group cross section table is then built only
        once and reused for every flux.

        It is recommended to make repeated calls to this method within a context
        manager using the :class:`openmc.lib.TemporarySession` class to avoid
        re-initializing OpenMC and loading cross sections each time.

        .. versionadded:: 0.15.0

        .. versionchanged:: 0.15.4
            ``multigroup_flux`` may be 2-D (or a list of 1-D arrays) to collapse
            several fluxes against a single shared cross section table, returning
            a list of :class:`MicroXS`. Added the ``cross_sections`` argument.

        Parameters
        ----------
        energies : iterable of float or str
            Energy group boundaries in [eV] or the name of the group structure
        multigroup_flux : iterable of float or iterable of iterable of float
            Energy-dependent multigroup flux values. Must be finite and
            non-negative. A 1-D input is a single flux; a 2-D input (or a list of
            1-D arrays) is a batch of fluxes that share the same group structure.
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
        cross_sections : PathLike, optional
            Cross section library used to resolve nuclide data availability and
            evaluate cross sections. Defaults to ``openmc.config['cross_sections']``.
        **init_kwargs : dict
            Keyword arguments passed to :func:`openmc.lib.init`

        Returns
        -------
        MicroXS or list of MicroXS
            A single :class:`MicroXS` for a 1-D ``multigroup_flux``; a list with
            one entry per flux for a 2-D input (a 1-row batch returns a
            1-element list, not an unwrapped :class:`MicroXS`).
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

        # A 1-D flux is a single domain; 2-D (or a list of 1-D) is a batch
        flux_array = np.asarray(multigroup_flux, dtype=float)
        if flux_array.ndim not in (1, 2):
            raise ValueError('multigroup_flux must be 1-D or 2-D')
        single = flux_array.ndim == 1
        fluxes = [flux_array] if single else list(flux_array)

        # check dimension consistency per flux
        n_groups = len(energies) - 1
        for flux in fluxes:
            if len(flux) != n_groups:
                raise ValueError('Length of flux array should be len(energies)-1')

        # Resolve the library once; data availability is derived from it
        if cross_sections is None:
            cross_sections = _find_cross_sections(model=None)
        nuclides_with_data = _get_nuclides_with_data(cross_sections)

        # Default nuclides/reactions from the chain only when needed
        if not nuclides or reactions is None:
            chain = _get_chain(chain_file)
            if not nuclides:
                nuclides = [nuc.name for nuc in chain.nuclides]
            if reactions is None:
                reactions = chain.reactions

        # Build the XS table once and collapse every flux against it
        table = _build_xs_table_ce(
            nuclides, reactions, energies, temperature, nuclides_with_data,
            cross_sections=cross_sections, **init_kwargs)
        micros = _collapse_fluxes(table, fluxes, nuclides, reactions)
        return micros[0] if single else micros

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
