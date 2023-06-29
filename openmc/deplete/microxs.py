"""MicroXS module

A pandas.DataFrame storing microscopic cross section data with
nuclide names as row indices and reaction names as column indices.
"""

import tempfile

from pandas import DataFrame, read_csv, Series
import numpy as np

from openmc.checkvalue import check_type, check_value, check_iterable_type
from openmc.exceptions import DataError
from openmc import StatePoint
import openmc
from .chain import Chain, REACTIONS
from .coupled_operator import _find_cross_sections, _get_nuclides_with_data

_valid_rxns = list(REACTIONS)
_valid_rxns.append('fission')


class MicroXS(DataFrame):
    """Microscopic cross section data for use in transport-independent depletion.

    .. versionadded:: 0.13.1

    """

    @classmethod
    def from_model(cls,
                   model,
                   domain,
                   nuclides=None,
                   reactions=None,
                   chain_file=None,
                   energy_bounds=(0, 20e6),
                   run_kwargs=None):
        """Generate a one-group cross-section dataframe using OpenMC.

        Note that the ``openmc`` executable must be compiled.

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
            Keyword arguments passed to :meth:`openmc.model.Model.run`

        Returns
        -------
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
        tallies = openmc.Tallies([rr_tally, flux_tally])

        model.tallies = tallies

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
        xs = reaction_rates / flux

        # Build Series objects
        series = {}
        for i, rx in enumerate(reactions):
            series[rx] = Series(xs[..., i], index=rr_tally.nuclides)

        # Revert to the original tallies and materials
        model.tallies = original_tallies

        return cls(series).rename_axis('nuclide')

    @classmethod
    def from_array(cls, nuclides, reactions, data):
        """
        Creates a ``MicroXS`` object from arrays.

        Parameters
        ----------
        nuclides : list of str
            List of nuclide symbols for that have data for at least one
            reaction.
        reactions : list of str
            List of reactions. All reactions must match those in
            :data:`openmc.deplete.chain.REACTIONS`
        data : ndarray of floats
            Array containing one-group microscopic cross section values for
            each nuclide and reaction. Cross section values are assumed to be
            in [b].

        Returns
        -------
        MicroXS
        """

        # Validate inputs
        if data.shape != (len(nuclides), len(reactions)):
            raise ValueError(
                f'Nuclides list of length {len(nuclides)} and '
                f'reactions array of length {len(reactions)} do not '
                f'match dimensions of data array of shape {data.shape}')

        cls._validate_micro_xs_inputs(
            nuclides, reactions, data)
        micro_xs = cls(index=nuclides, columns=reactions, data=data)

        return micro_xs

    @classmethod
    def from_csv(cls, csv_file, **kwargs):
        """
        Load a ``MicroXS`` object from a ``.csv`` file.

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

        micro_xs = cls(read_csv(csv_file, index_col=0, **kwargs))

        cls._validate_micro_xs_inputs(list(micro_xs.index),
                                      list(micro_xs.columns),
                                      micro_xs.to_numpy())
        return micro_xs

    @staticmethod
    def _validate_micro_xs_inputs(nuclides, reactions, data):
        check_iterable_type('nuclides', nuclides, str)
        check_iterable_type('reactions', reactions, str)
        check_type('data', data, np.ndarray, expected_iter_type=float)
        for reaction in reactions:
            check_value('reactions', reaction, _valid_rxns)
