"""MicroXS module

A pandas.DataFrame storing microscopic cross section data with
nuclides names as row indices and reaction names as column indices.
"""

import tempfile
from pathlib import Path

from pandas import DataFrame, read_csv, concat
from numpy import ndarray

from openmc.checkvalue import check_type, check_value, check_iterable_type
from openmc.mgxs import EnergyGroups, ArbitraryXS, FissionXS
from openmc import Tallies, StatePoint

from .chain import REACTIONS

_valid_rxns = list(REACTIONS)
_valid_rxns.append('fission')


class MicroXS(DataFrame):
    """Stores microscopic cross section data for use in
    independent depletion.
    """

    @classmethod
    def from_model(cls,
                   model,
                   reaction_domain,
                   reactions=['(n,gamma)',
                              '(n,2n)',
                              '(n,p)',
                              '(n,a)',
                              '(n,3n)',
                              '(n,4n)',
                              'fission'],
                   energy_bounds=(0, 20e6)):
        """Generate a one-group cross-section dataframe using
        OpenMC. Note that the ``openmc`` executable must be compiled.

        Parameters
        ----------
        model : openmc.Model
            OpenMC model object. Must contain geometry, materials, and settings.
        reaction_domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
            Domain in which to tally reaction rates.
        reactions : list of str, optional
            Reaction names to tally
        energy_bound : 2-tuple of float, optional
            Bounds for the energy group.

        Returns
        -------
        MicroXS
            Cross section data in [b]

        """
        groups = EnergyGroups(energy_bounds)

        # Set up the reaction tallies
        original_tallies = model.tallies
        tallies = Tallies()
        xs = {}
        for rxn in reactions:
            if rxn == 'fission':
                xs[rxn] = FissionXS(domain=reaction_domain, groups=groups, by_nuclide=True)
            else:
                xs[rxn] = ArbitraryXS(rxn, domain=reaction_domain, groups=groups, by_nuclide=True)
            tallies += xs[rxn].tallies.values()

        model.tallies = tallies

        # create temporary run
        with tempfile.TemporaryDirectory() as temp_dir:
            statepoint_path = model.run(cwd=temp_dir)

            with StatePoint(statepoint_path) as sp:
                for rxn in xs:
                    xs[rxn].load_from_statepoint(sp)

        # Build the DataFrame
        micro_xs = cls()
        for rxn in xs:
            df = xs[rxn].get_pandas_dataframe(xs_type='micro')
            df.index = df['nuclide']
            df.drop(['nuclide', xs[rxn].domain_type, 'group in', 'std. dev.'], axis=1, inplace=True)
            df.rename({'mean': rxn}, axis=1, inplace=True)
            micro_xs = concat([micro_xs, df], axis=1)

        micro_xs._units = 'b'

        # Revert to the original tallies
        model.tallies = original_tallies

        return micro_xs

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
        check_type('data', data, ndarray, expected_iter_type=float)
        for reaction in reactions:
            check_value('reactions', reaction, _valid_rxns)
