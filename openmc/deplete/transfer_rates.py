from numbers import Real
import re

from openmc.checkvalue import check_type, check_value
from openmc import Material
from openmc.data import ELEMENT_SYMBOL


class TransferRates:
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
    model : openmc.Model
        OpenMC model containing materials and geometry. If using
        :class:`openmc.deplete.CoupledOperator`, the model must also contain
        a :class:`opnemc.Settings` object.

    Attributes
    ----------
    burnable_mats : list of str
        All burnable material IDs.
    local_mats : list of str
        All burnable material IDs being managed by a single process
    transfer_rates : dict of str to dict
        Container of transfer rates, components (elements and/or nuclides) and
        destination material
    index_transfer : Set of pair of str
        Pair of strings needed to build final matrix (destination_material, mat)
    """

    def __init__(self, operator, model):

        self.materials = model.materials
        self.burnable_mats = operator.burnable_mats
        self.local_mats = operator.local_mats

        # initialize transfer rates container dict
        self.transfer_rates = {mat: {} for mat in self.burnable_mats}
        self.index_transfer = set()

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
            check_value("Depeletable Material", str(val.id), self.burnable_mats)
            val = val.id

        elif isinstance(val, str):
            if val.isnumeric():
                check_value("Material ID", str(val), self.burnable_mats)
            else:
                check_value(
                    "Material name",
                    val,
                    [mat.name for mat in self.materials if mat.depletable],
                )
                val = [mat.id for mat in self.materials if mat.name == val][0]

        elif isinstance(val, int):
            check_value("Material ID", str(val), self.burnable_mats)

        return str(val)

    def get_transfer_rate(self, material, component):
        """Return transfer rate for given material and element.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        component : str
            Element or nuclide to get transfer rate value

        Returns
        -------
        transfer_rate : list of floats
            Transfer rate values

        """
        material_id = self._get_material_id(material)
        check_type("component", component, str)
        return [i[0] for i in self.transfer_rates[material_id][component]]

    def get_destination_material(self, material, component):
        """Return destination material for given material and
        component, if defined.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        component : str
            Element or nuclide that gets transferred to another material.

        Returns
        -------
        destination_material_id : list of str
            Depletable material ID to where the element or nuclide gets
            transferred

        """
        material_id = self._get_material_id(material)
        check_type("component", component, str)
        if component in self.transfer_rates[material_id]:
            return [i[1] for i in self.transfer_rates[material_id][component]]
        else:
            return []

    def get_components(self, material):
        """Extract removing elements and/or nuclides for a given material

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material

        Returns
        -------
        elements : list
            List of elements and nuclides where transfer rates exist

        """
        material_id = self._get_material_id(material)
        if material_id in self.transfer_rates:
            return self.transfer_rates[material_id].keys()

    def set_transfer_rate(
        self,
        material,
        components,
        transfer_rate,
        transfer_rate_units="1/s",
        destination_material=None,
    ):
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
        destination_material : openmc.Material or str or int, Optional
            Destination material to where nuclides get fed.
        transfer_rate_units : {'1/s', '1/min', '1/h', '1/d', '1/a'}
            Units for values specified in the transfer_rate argument. 's' for
            seconds, 'min' for minutes, 'h' for hours, 'a' for Julian years.

        """
        material_id = self._get_material_id(material)
        check_type("transfer_rate", transfer_rate, Real)
        check_type("components", components, list, expected_iter_type=str)

        if destination_material is not None:
            destination_material_id = self._get_material_id(destination_material)
            if len(self.burnable_mats) > 1:
                check_value(
                    "destination_material",
                    str(destination_material_id),
                    self.burnable_mats,
                )
            else:
                raise ValueError(
                    "Transfer to material "
                    f"{destination_material_id} is set, but there "
                    "is only one depletable material"
                )
        else:
            destination_material_id = None

        if transfer_rate_units in ("1/s", "1/sec"):
            unit_conv = 1
        elif transfer_rate_units in ("1/min", "1/minute"):
            unit_conv = 60
        elif transfer_rate_units in ("1/h", "1/hr", "1/hour"):
            unit_conv = 60 * 60
        elif transfer_rate_units in ("1/d", "1/day"):
            unit_conv = 24 * 60 * 60
        elif transfer_rate_units in ("1/a", "1/year"):
            unit_conv = 365.25 * 24 * 60 * 60
        else:
            raise ValueError("Invalid transfer rate unit " f'"{transfer_rate_units}"')

        for component in components:
            current_components = self.transfer_rates[material_id].keys()
            split_component = re.split(r"\d+", component)
            element = split_component[0]
            if element not in ELEMENT_SYMBOL.values():
                raise ValueError(f"{component} is not a valid nuclide or " "element.")
            else:
                if len(split_component) == 1:
                    element_nucs = [
                        c for c in current_components if re.match(component + r"\d", c)
                    ]
                    if len(element_nucs) > 0:
                        nuc_str = ", ".join(element_nucs)
                        raise ValueError(
                            "Cannot add transfer rate for element "
                            f"{component} to material {material_id} "
                            f"with transfer rate(s) for nuclide(s) "
                            f"{nuc_str}."
                        )

                else:
                    if element in current_components:
                        raise ValueError(
                            "Cannot add transfer rate for nuclide "
                            f"{component} to material {material_id} "
                            f"where element {element} already has "
                            "a transfer rate."
                        )

            if component in self.transfer_rates[material_id]:
                self.transfer_rates[material_id][component].append(
                    (transfer_rate / unit_conv, destination_material_id)
                )
            else:
                self.transfer_rates[material_id][component] = [
                    (transfer_rate / unit_conv, destination_material_id)
                ]
            if destination_material_id is not None:
                self.index_transfer.add((destination_material_id, material_id))
