"""MicroXS module

A class for storing microscopic cross section data that can be used with the
IndependentOperator class for depletion.
"""

from __future__ import annotations
from collections.abc import Sequence
import shutil
from tempfile import TemporaryDirectory
from typing import Union, TypeAlias

import pandas as pd
import numpy as np

from openmc.checkvalue import check_type, check_value, check_iterable_type, PathLike
from openmc.utility_funcs import change_directory
from openmc import StatePoint
from openmc.mgxs import GROUP_STRUCTURES
from openmc.data import REACTION_MT
import openmc
from .chain import Chain, REACTIONS, _get_chain
from .coupled_operator import _find_cross_sections, _get_nuclides_with_data
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
    run_kwargs=None
) -> tuple[list[np.ndarray], list[MicroXS]]:
    """Generate microscopic cross sections and fluxes for multiple domains.

    This function runs a neutron transport solve to obtain the flux and reaction
    rates in the specified domains and computes multigroup microscopic cross
    sections that can be used in depletion calculations with the
    :class:`~openmc.deplete.IndependentOperator` class.

    .. versionadded:: 0.14.0

    .. versionchanged:: 0.15.3
        Added `reaction_rate_mode` and `path_statepoint` arguments.

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
        Indicate how reaction rates should be calculated. The "direct" method
        tallies reaction rates directly. The "flux" method tallies a multigroup
        flux spectrum and then collapses multigroup reaction rates after a
        transport solve (with an option to tally some reaction rates directly).
    chain_file : PathLike or Chain, optional
        Path to the depletion chain XML file or an instance of
        openmc.deplete.Chain. Used to determine cross sections for materials not
        present in the inital composition. Defaults to
        ``openmc.config['chain_file']``.
    path_statepoint : path-like, optional
        Path to write the statepoint file from the neutron transport solve to.
        By default, The statepoint file is written to a temporary directory and
        is not kept.
    run_kwargs : dict, optional
        Keyword arguments passed to :meth:`openmc.Model.run`

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
    original_tallies = model.tallies

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

    if reaction_rate_mode == 'direct':
        rr_tally = openmc.Tally(name='MicroXS RR')
        rr_tally.filters = [domain_filter, energy_filter]
        rr_tally.nuclides = nuclides
        rr_tally.multiply_density = False
        rr_tally.scores = reactions
        model.tallies.append(rr_tally)

    if openmc.lib.is_initialized:
        openmc.lib.finalize()

        if comm.rank == 0:
            model.export_to_model_xml()
        comm.barrier()
        # Reinitialize with tallies
        openmc.lib.init(intracomm=comm)

    # create temporary run
    with TemporaryDirectory() as temp_dir:
        if run_kwargs is None:
            run_kwargs = {}
        else:
            run_kwargs = dict(run_kwargs)
        run_kwargs.setdefault('cwd', temp_dir)
        statepoint_path = model.run(**run_kwargs)

        if comm.rank == 0:
            # Move the statepoint file if it is being saved to a specific path
            if path_statepoint is not None:
                shutil.move(statepoint_path, path_statepoint)
                statepoint_path = path_statepoint

            with StatePoint(statepoint_path) as sp:
                if reaction_rate_mode == 'direct':
                    rr_tally = sp.tallies[rr_tally.id]
                    rr_tally._read_results()
                flux_tally = sp.tallies[flux_tally.id]
                flux_tally._read_results()

    # Get flux values and make energy groups last dimension
    flux_tally = comm.bcast(flux_tally)
    flux = flux_tally.get_reshaped_data()  # (domains, groups, 1, 1)
    flux = np.moveaxis(flux, 1, -1)  # (domains, 1, 1, groups)

    # Create list where each item corresponds to one domain
    fluxes = list(flux.squeeze((1, 2)))

    if reaction_rate_mode == 'direct':
        # Get reaction rates
        rr_tally = comm.bcast(rr_tally)
        reaction_rates = rr_tally.get_reshaped_data()  # (domains, groups, nuclides, reactions)

        # Make energy groups last dimension
        reaction_rates = np.moveaxis(reaction_rates, 1, -1)  # (domains, nuclides, reactions, groups)

        # Divide RR by flux to get microscopic cross sections. The indexing
        # ensures that only non-zero flux values are used, and broadcasting is
        # applied to align the shapes of reaction_rates and flux for division.
        xs = np.empty_like(reaction_rates) # (domains, nuclides, reactions, groups)
        d, _, _, g = np.nonzero(flux)
        xs[d, ..., g] = reaction_rates[d, ..., g] / flux[d, :, :, g]

        # Create lists where each item corresponds to one domain
        micros = [MicroXS(xs_i, nuclides, reactions) for xs_i in xs]
    else:
        micros = [MicroXS.from_multigroup_flux(
            energies=energies,
            multigroup_flux=flux_i,
            chain_file=chain_file,
            nuclides=nuclides,
            reactions=reactions
        ) for flux_i in fluxes]

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

        # Normalize multigroup flux
        multigroup_flux = np.array(multigroup_flux)
        multigroup_flux /= multigroup_flux.sum()

        # Create 3D array for microscopic cross sections
        microxs_arr = np.zeros((len(nuclides), len(mts), 1))

        # Create a material with all nuclides
        mat_all_nucs = openmc.Material()
        for nuc in nuclides:
            if nuc in nuclides_with_data:
                mat_all_nucs.add_nuclide(nuc, 1.0)
        mat_all_nucs.set_density("atom/b-cm", 1.0)

        # Create simple model containing the above material
        surf1 = openmc.Sphere(boundary_type="vacuum")
        surf1_cell = openmc.Cell(fill=mat_all_nucs, region=-surf1)
        model = openmc.Model()
        model.geometry = openmc.Geometry([surf1_cell])
        model.settings = openmc.Settings(
            particles=1, batches=1, output={'summary': False})

        with change_directory(tmpdir=True):
            # Export model within temporary directory
            model.export_to_model_xml()

            with openmc.lib.run_in_memory(**init_kwargs):
                # For each nuclide and reaction, compute the flux-averaged
                # cross section
                for nuc_index, nuc in enumerate(nuclides):
                    if nuc not in nuclides_with_data:
                        continue
                    lib_nuc = openmc.lib.nuclides[nuc]
                    for mt_index, mt in enumerate(mts):
                        xs = lib_nuc.collapse_rate(
                            mt, temperature, energies, multigroup_flux
                        )
                        microxs_arr[nuc_index, mt_index, 0] = xs

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
