from numbers import Real
import re

from openmc.checkvalue import check_type, check_value
from openmc import Material
from openmc.data import ELEMENT_SYMBOL, isotopes, AVOGADRO, atomic_mass

from math import isclose
import numpy as np

class ExternalRates:
    """External rates class for defining addition terms of depletion equation.

    .. versionadded:: 0.15.1

    Parameters
    ----------
    operator : openmc.TransportOperator
        Depletion operator
    materials : openmc.Model
        OpenMC materials.
    number_of_timesteps : int
        Total number of depletion timesteps

    Attributes
    ----------
    burnable_mats : list of str
        All burnable material IDs.
    local_mats : list of str
        All burnable material IDs being managed by a single process
    number_of_timesteps : int
        Total number of depletion timesteps
    external_rates : dict of str to dict
        Container of timesteps, external rates, components (elements and/or
        nuclides) and optionally destination material
    external_timesteps : list of int
        Container of all timesteps indeces with an external rate defined.
    """

    def __init__(self, operator, materials, number_of_timesteps):

        self.materials = materials
        self.burnable_mats = operator.burnable_mats
        self.local_mats = operator.local_mats
        self.number_of_timesteps = number_of_timesteps

        #initialize transfer rates container dict
        self.external_rates = {mat: {} for mat in self.burnable_mats}
        self.external_timesteps =  []

    def _get_material_id(self, val):
        """Helper method for getting material id from Material obj or name.

        Parameters
        ----------
        val : openmc.Material or str or int representing material name/id

        Returns
        -------
        material_id : str

        """
        if isinstance(val, Material):
            check_value('Depeletable Material', str(val.id), self.burnable_mats)
            val = val.id

        elif isinstance(val, str):
            if val.isnumeric():
                check_value('Material ID', str(val), self.burnable_mats)
            else:
                check_value('Material name', val,
                        [mat.name for mat in self.materials if mat.depletable])
                val = [mat.id for mat in self.materials if mat.name == val][0]

        elif isinstance(val, int):
            check_value('Material ID', str(val), self.burnable_mats)

        return str(val)

    def get_external_rate(self, material, component, timestep,
                          destination_material=None):
        """Return transfer rate for given material and element.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        component : str
            Element or nuclide to get transfer rate value
        timestep : int
            Current timestep index
        destination_material : openmc.Material or str or int, Optional
            Destination material to where nuclides get fed

        Returns
        -------
        external_rate : list of floats
            External rate values

        """
        material_id = self._get_material_id(material)
        check_type('component', component, str)
        if destination_material is not None:
            dest_mat_id = self._get_material_id(destination_material)
            return [i[1] for i in self.external_rates[material_id][component]
                if timestep in i[0] and dest_mat_id==i[2]]
        else:
            return [i[1] for i in self.external_rates[material_id][component]
                if timestep in i[0]]

    def get_components(self, material, timestep, destination_material=None):
        """Extract removing elements and/or nuclides for a given material at a
        given timestep

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        timestep : int
            Current timestep index
        destination_material : openmc.Material or str or int, Optional
            Destination material to where nuclides get fed

        Returns
        -------
        components : list
            List of elements or nuclides with external rates set at a given
            timestep

        """
        material_id = self._get_material_id(material)
        if destination_material is not None:
            dest_mat_id = self._get_material_id(destination_material)
        else:
            dest_mat_id = None

        all_components = []
        if material_id in self.external_rates:
            mat_components = self.external_rates[material_id]

            for component in mat_components:
                if dest_mat_id:
                    #check for both timestep and destination material ids
                    if np.isin(timestep,
                            [val[0] for val in mat_components[component]]) and \
                       np.isin(dest_mat_id,
                            [val[2] for val in mat_components[component]]):
                        all_components.append(component)
                else:
                    #check only for timesteps
                    if np.isin(timestep,
                            [val[0] for val in mat_components[component]]):
                        all_components.append(component)
        return all_components

class TransferRates(ExternalRates):
    """Class for defining continuous removals and feeds.

    Molten Salt Reactors (MSRs) benefit from continuous reprocessing,
    which removes fission products and feeds fresh fuel into the system. MSRs
    inspired the development of this class.

    An instance of this class can be passed directly to an instance of one of
    the :class:`openmc.deplete.Integrator` classes.

    .. versionadded:: 0.14.0

    Parameters
    ----------
    operator : openmc.TransportOperator
        Depletion operator
    materials : openmc.Model
        OpenMC materials.
    number_of_timesteps : int
        Total number of depletion timesteps

    Attributes
    ----------
    burnable_mats : list of str
        All burnable material IDs.
    local_mats : list of str
        All burnable material IDs being managed by a single process
    external_rates : dict of str to dict
        Container of timesteps, transfer rates, components (elements and/or
        nuclides) and destination material
    external_timesteps : list of int
        Container of all timesteps indeces with an external rate defined.
    index_transfer : Set of pair of str
        Pair of strings needed to build final matrix (destination_material, mat)
    """

    def __init__(self, operator, materials, number_of_timesteps):

        super().__init__(
            operator,
            materials,
            number_of_timesteps
        )

        self.index_transfer = dict()
        self.chain_nuclides = [nuc.name for nuc in operator.chain.nuclides]

    def set_transfer_rate(self, material, components, transfer_rate,
                          transfer_rate_units='1/s', timesteps=None,
                          destination_material=None):
        """Set element and/or nuclide transfer rates in a depletable material.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        components : list of str
            List of strings of elements and/or nuclides that share transfer rate.
            Cannot add transfer rates for nuclides to a material where a
            transfer rate for its element is specified and vice versa.
        transfer_rate : float
            Rate at which elements and/or nuclides are transferred. A positive or
            negative value corresponds to a removal or feed rate, respectively.
        transfer_rate_units : {'1/s', '1/min', '1/h', '1/d', '1/a'}
            Units for values specified in the transfer_rate argument. 's' for
            seconds, 'min' for minutes, 'h' for hours, 'a' for Julian years.
        timesteps : list of int, Optional
            List of timestep indeces where to set transfer rates.
            Default to None means the transfer rate is set for all timesteps.
        destination_material : openmc.Material or str or int, Optional
            Destination material to where nuclides get fed.

        """
        material_id = self._get_material_id(material)
        check_type('transfer_rate', transfer_rate, Real)
        check_type('components', components, list, expected_iter_type=str)

        if destination_material is not None:
            destination_material_id = self._get_material_id(destination_material)
            if len(self.burnable_mats) > 1:
                check_value('destination_material', str(destination_material_id),
                            self.burnable_mats)
            else:
                raise ValueError('Transfer to material '
                                 f'{destination_material_id} is set, but there '
                                 'is only one depletable material')
        else:
            destination_material_id = None

        if transfer_rate_units in ('1/s', '1/sec'):
            unit_conv = 1
        elif transfer_rate_units in ('1/min', '1/minute'):
            unit_conv = 60
        elif transfer_rate_units in ('1/h', '1/hr', '1/hour'):
            unit_conv = 60*60
        elif transfer_rate_units in ('1/d', '1/day'):
            unit_conv = 24*60*60
        elif transfer_rate_units in ('1/a', '1/year'):
            unit_conv = 365.25*24*60*60
        else:
            raise ValueError('Invalid transfer rate unit '
                             f'"{transfer_rate_units}"')

        if timesteps is not None:
            for timestep in timesteps:
                check_value('timestep', timestep, range(self.number_of_timesteps))
            timesteps = np.array(timesteps)
        else:
            timesteps = np.arange(self.number_of_timesteps)

        for component in components:
            current_components = self.external_rates[material_id].keys()
            split_component = re.split(r'\d+', component)
            element = split_component[0]
            if element not in ELEMENT_SYMBOL.values():
                raise ValueError(f'{component} is not a valid nuclide or '
                                 'element.')
            else:
                if len(split_component) == 1:
                    element_nucs = [c for c in current_components
                                    if re.match(component + r'\d', c)]
                    if len(element_nucs) > 0:
                        nuc_str = ", ".join(element_nucs)
                        raise ValueError('Cannot add transfer rate for element '
                                         f'{component} to material {material_id} '
                                         f'with transfer rate(s) for nuclide(s) '
                                         f'{nuc_str}.')

                else:
                    if element in current_components:
                        raise ValueError('Cannot add transfer rate for nuclide '
                                         f'{component} to material {material_id} '
                                         f'where element {element} already has '
                                         'a transfer rate.')

            if component in self.external_rates[material_id]:
                self.external_rates[material_id][component].append(
                    (timesteps,
                     transfer_rate / unit_conv,
                     destination_material_id))
            else:
                self.external_rates[material_id][component] = [
                    (timesteps,
                     transfer_rate / unit_conv,
                     destination_material_id)]

            if destination_material_id is not None:
                for timestep in timesteps:
                    if timestep not in self.index_transfer:
                        self.index_transfer[timestep] = [(destination_material_id,
                                                        material_id)]
                    else:
                        self.index_transfer[timestep].append((destination_material_id,
                                                        material_id))

            self.external_timesteps = np.unique(np.concatenate(
                    [self.external_timesteps, timesteps]))


class ExternalSourceRates(ExternalRates):
    """Class for defining external source rates.

    An instance of this class can be passed directly to an instance of one of
    the :class:`openmc.deplete.Integrator` classes.

    .. versionadded:: 0.15.1

    Parameters
    ----------
    operator : openmc.TransportOperator
        Depletion operator
    materials : openmc.Model
        OpenMC materials.
    number_of_timesteps : int
        Total number of depletion timesteps

    Attributes
    ----------
    burnable_mats : list of str
        All burnable material IDs.
    local_mats : list of str
        All burnable material IDs being managed by a single process
    external_timesteps : list of int
        Container of all timesteps indeces with an external rate defined.
    external_rates : dict of str to dict
        Container of timesteps external source rates, and components
        (elements and/or nuclides)
    """

    def __init__(self, operator, materials, number_of_timesteps):

        super().__init__(
            operator,
            materials,
            number_of_timesteps
        )

    def reformat_nuclide_vectors(self, vectors):
        """Remove last element of nuclide vector that was added for homogenize
        the vector by the depletion solver.

        Parameters
        ----------
        vectors : list of array
            List of nuclides vector to reformat

        Returns
        -------
        vectors : list of array
            Updated list of nuclides vector

        """
        mat_indeces = [idx for idx,i in enumerate(self.local_mats) \
                        if self.external_rates[i]]

        for mat_index in mat_indeces:
            vectors[mat_index] = vectors[mat_index][:-1]

        return vectors

    def set_external_source_rate(self, material, external_source_vector,
            external_source_rate, external_source_rate_units='g/s', timesteps=None):
        """Set element and/or nuclide composition vector external source rates
        to a depletable material.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        external_source_vector : dict of str to float
            External source rate composition vector, where key can be an element
            or a nuclide and value the corresponding weigth percent.
        external_source_rate : float
            External source rate in unit of grams per time. A positive or
            negative value corresponds to a feed or removal rate, respectively.
        external_source_rate_units : {'g/s', 'g/min', 'g/h', 'g/d', 'g/a'}
            Units for values specified in the external_source_rate argument.
            's' for seconds, 'min' for minutes, 'h' for hours, 'a' for
            Julian years.
        timesteps : list of int, Optional
            List of timestep indeces where to set external source rates.
            Default to None means the external source rate is set for all
            timesteps.

        """

        material_id = self._get_material_id(material)
        check_type('external_source_rate', external_source_rate, Real)
        check_type('external_source_vector', external_source_vector, dict,
                        expected_iter_type=str)

        if external_source_rate_units in ('g/s', 'g/sec'):
            unit_conv = 1
        elif external_source_rate_units in ('g/min', 'g/minute'):
            unit_conv = 60
        elif external_source_rate_units in ('g/h', 'g/hr', 'g/hour'):
            unit_conv = 60*60
        elif external_source_rate_units in ('g/d', 'g/day'):
            unit_conv = 24*60*60
        elif external_source_rate_units in ('g/a', 'g/year'):
            unit_conv = 365.25*24*60*60
        else:
            raise ValueError('Invalid external source rate unit '
                             f'"{external_source_rate_units}"')

        if timesteps is not None:
            for timestep in timesteps:
                check_value('timestep', timestep, range(self.number_of_timesteps))
            timesteps = np.array(timesteps)
        else:
            timesteps = np.arange(self.number_of_timesteps)


        percents = external_source_vector.values()
        norm_percents = [float(i) / sum(percents) for i in percents]
        components = external_source_vector.keys()

        atoms_per_nuc = {}
        for component, percent in zip(components, norm_percents):
            split_component = re.split(r'\d+', component)
            element = split_component[0]
            if element not in ELEMENT_SYMBOL.values():
                raise ValueError(f'{component} is not a valid nuclide or '
                                 'element.')
            else:
                if len(split_component) == 1:
                    if not isotopes(component):
                        raise ValueError(f'Cannot add element {component} '
                                         'as it is not naturally abundant. '
                                         'Specify a nuclide vector instead. ')
                    else:
                        for nuc, frac in isotopes(component):
                            atoms_per_nuc[nuc] = \
                                external_source_rate / atomic_mass(nuc) \
                                * AVOGADRO * frac * percent \
                                / unit_conv

                else:
                    atoms_per_nuc[component] = \
                            external_source_rate / atomic_mass(component) \
                            * AVOGADRO * percent / unit_conv

        for nuc, val in atoms_per_nuc.items():
            if nuc in self.external_rates[material_id]:
                self.external_rates[material_id][nuc].append((timesteps, val, None))
            else:
                self.external_rates[material_id][nuc] = [(timesteps, val, None)]

        self.external_timesteps = np.unique(np.concatenate(
                    [self.external_timesteps, timesteps]))
