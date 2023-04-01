"""MicroXS module

A pandas.DataFrame storing microscopic cross section data with
nuclide names as row indices and reaction names as column indices.
"""

import tempfile
from copy import deepcopy

from pandas import DataFrame, read_csv
import numpy as np

from openmc.checkvalue import check_type, check_value, check_iterable_type
from openmc.exceptions import DataError
from openmc.mgxs import EnergyGroups, ArbitraryXS, FissionXS
from openmc import Tallies, StatePoint, Materials
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
                   reaction_domain,
                   chain_file=None,
                   dilute_initial=1.0e3,
                   energy_bounds=(0, 20e6),
                   run_kwargs=None):
        """Generate a one-group cross-section dataframe using
        OpenMC. Note that the ``openmc`` executable must be compiled.

        Parameters
        ----------
        model : openmc.Model
            OpenMC model object. Must contain geometry, materials, and settings.
        reaction_domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
            Domain in which to tally reaction rates.
        chain_file : str, optional
            Path to the depletion chain XML file that will be used in depletion
            simulation. Used to determine cross sections for materials not
            present in the inital composition. Defaults to
            ``openmc.config['chain_file']``.
        dilute_initial : float, optional
            Initial atom density [atoms/cm^3] to add for nuclides that
            are zero in initial condition to ensure they exist in the cross
            section data. Only done for nuclides with reaction rates.
        reactions : list of str, optional
            Reaction names to tally
        energy_bound : 2-tuple of float, optional
            Bounds for the energy group.
        run_kwargs : dict, optional
            Keyword arguments passed to :meth:`openmc.model.Model.run`

        Returns
        -------
        MicroXS
            Cross section data in [b]

        """
        groups = EnergyGroups(energy_bounds)

        # Set up the reaction tallies
        original_tallies = model.tallies
        original_materials = deepcopy(model.materials)
        tallies = Tallies()
        xs = {}
        reactions, diluted_materials = cls._add_dilute_nuclides(chain_file,
                                                                model,
                                                                dilute_initial)
        model.materials = diluted_materials

        for rx in reactions:
            if rx == 'fission':
                xs[rx] = FissionXS(domain=reaction_domain,
                                   energy_groups=groups, by_nuclide=True)
            else:
                xs[rx] = ArbitraryXS(rx, domain=reaction_domain,
                                     energy_groups=groups, by_nuclide=True)
            tallies += xs[rx].tallies.values()

        model.tallies = tallies

        # create temporary run
        with tempfile.TemporaryDirectory() as temp_dir:
            if run_kwargs is None:
                run_kwargs = {}
            run_kwargs.setdefault('cwd', temp_dir)
            statepoint_path = model.run(**run_kwargs)

            with StatePoint(statepoint_path) as sp:
                for rx in xs:
                    xs[rx].load_from_statepoint(sp)

        # Build the DataFrame
        series = {}
        for rx in xs:
            df = xs[rx].get_pandas_dataframe(xs_type='micro')
            series[rx] = df.set_index('nuclide')['mean']

        # Revert to the original tallies and materials
        model.tallies = original_tallies
        model.materials = original_materials

        return cls(series)

    @classmethod
    def _add_dilute_nuclides(cls, chain_file, model, dilute_initial):
        """
        Add nuclides not present in burnable materials that have neutron data
        and are present in the depletion chain to those materials. This allows
        us to tally those specific nuclides for reactions to create one-group
        cross sections.

        Parameters
        ----------
        chain_file : str
            Path to the depletion chain XML file that will be used in depletion
            simulation. Used to determine cross sections for materials not
            present in the inital composition.
        model : openmc.Model
            Model object
        dilute_initial : float
            Initial atom density [atoms/cm^3] to add for nuclides that
            are zero in initial condition to ensure they exist in the cross
            section data. Only done for nuclides with reaction rates.

        Returns
        -------
        reactions : list of str
            List of reaction names
        diluted_materials : openmc.Materials
            :class:`openmc.Materials` object with nuclides added to burnable
            materials.
        """
        if chain_file is None:
            chain_file = openmc.config.get('chain_file')
            if chain_file is None:
                raise DataError(
                    "No depletion chain specified and could not find depletion "
                    "chain in openmc.config['chain_file']"
                )
        chain = Chain.from_xml(chain_file)
        reactions = chain.reactions
        cross_sections = _find_cross_sections(model)
        nuclides_with_data = _get_nuclides_with_data(cross_sections)
        burnable_nucs = [nuc.name for nuc in chain.nuclides
                         if nuc.name in nuclides_with_data]
        diluted_materials = Materials()
        for material in model.materials:
            if material.depletable:
                nuc_densities = material.get_nuclide_atom_densities()
                dilute_density = 1.0e-24 * dilute_initial
                material.set_density('sum')
                for nuc, density in nuc_densities.items():
                    material.remove_nuclide(nuc)
                    material.add_nuclide(nuc, density)
                for burn_nuc in burnable_nucs:
                    if burn_nuc not in nuc_densities:
                        material.add_nuclide(burn_nuc,
                                                     dilute_density)
            diluted_materials.append(material)

        return reactions, diluted_materials

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
