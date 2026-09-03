import lxml.etree as ET

import openmc
import openmc.checkvalue as cv
from .field import VelocityField


class DNPDrift:
    """Settings for delayed neutron precursor drift.

    This class defines the parameters needed to model the transport of
    delayed neutron precursors in a flowing fuel system (e.g., molten
    salt reactor).

    Parameters
    ----------
    velocity_field : openmc.VelocityField
        Velocity field describing the fluid motion.
    boundary_map : dict
        Mapping of boundary types to physical groups.
        Keys: "inlet", "outlet", "wall" (optional). Values: list of int.
        "inlet" and "outlet" are required.
    mapping : {'nodal', 'cell'}
        How field values map to the mesh.
    integrator : {'RK4'}
        Time integration scheme for precursor transport.
    integrator_dt : float
        Time step for integration in seconds.
    recycling : bool
        Whether precursors are recycled from outlet to inlet.
    external_travel_time : float
        Time for precursors to travel through the external loop in
        seconds. Used when recycling is True. Default is 1.0s.

    Attributes
    ----------
    velocity_field : openmc.VelocityField
        Velocity field describing the fluid motion.
    boundary_map : dict
        Mapping of boundary types to physical groups.
    mapping : {'nodal', 'cell'}
        How field values map to the mesh.
    integrator : {'RK4'}
        Time integration scheme.
    integrator_dt : float
        Time step for integration in seconds.
    recycling : bool
        Whether precursors are recycled.
    external_travel_time : float
        External travel time in seconds.

    """
    def __init__(
        self,
        velocity_field,
        boundary_map,
        mapping="cell",
        integrator="RK4",
        integrator_dt=0.1,
        recycling=True,
        external_travel_time=1.0,
    ):
        self.velocity_field = velocity_field
        self.boundary_map = boundary_map
        self.mapping = mapping
        self.integrator = integrator
        self.integrator_dt = integrator_dt
        self.recycling = recycling
        self.external_travel_time = external_travel_time

    @property
    def velocity_field(self) -> VelocityField:
        return self._velocity_field

    @velocity_field.setter
    def velocity_field(self, velocity_field: VelocityField):
        cv.check_type("velocity_field", velocity_field, VelocityField)
        self._velocity_field = velocity_field

    @property
    def boundary_map(self) -> dict[str, list[int]]:
        return self._boundary_map

    @boundary_map.setter
    def boundary_map(self, boundary_map: dict[str, int | list[int]]):
        cv.check_type("boundary_map", boundary_map, dict)

        # Create a copy to avoid mutating the original object
        bm = {k: v if isinstance(v, list) else [v] for k, v in boundary_map.items()}

        # Check that inlet and outlet are defined
        required_keys = {"inlet", "outlet"}
        missing_keys = required_keys - bm.keys()
        if missing_keys:
            raise ValueError(
                f"boundary_map is missing required keys: {missing_keys}"
            )

        # Check entries
        valid_keys = {"inlet", "outlet", "wall"}
        for key in bm:
            if key not in valid_keys:
                raise ValueError(
                    f"Invalid boundary type '{key}'. "
                    f"Must be one of {valid_keys}."
                )
            cv.check_type(f"boundary_map[{key}]", bm[key], list)
            for sid in bm[key]:
                cv.check_type(f"surface id in boundary_map[{key}]", sid, int)
        self._boundary_map = bm

    @property
    def mapping(self) -> str:
        return self._mapping

    @mapping.setter
    def mapping(self, mapping: str):
        cv.check_value("mapping", mapping, ("nodal", "cell"))
        self._mapping = mapping

    @property
    def integrator(self) -> str:
        return self._integrator

    @integrator.setter
    def integrator(self, integrator: str):
        cv.check_value("integrator", integrator, ("RK4",))
        self._integrator = integrator

    @property
    def integrator_dt(self) -> float:
        return self._integrator_dt

    @integrator_dt.setter
    def integrator_dt(self, integrator_dt: int | float):
        cv.check_type("integrator_dt", integrator_dt, (int, float))
        cv.check_greater_than("integrator_dt", integrator_dt, 0)
        self._integrator_dt = float(integrator_dt)

    @property
    def recycling(self) -> bool:
        return self._recycling

    @recycling.setter
    def recycling(self, recycling: bool):
        cv.check_type("recycling", recycling, bool)
        self._recycling = recycling

    @property
    def external_travel_time(self) -> float:
        return self._external_travel_time

    @external_travel_time.setter
    def external_travel_time(self, ett: int | float):
        cv.check_type("external_travel_time", ett, (int, float))
        cv.check_greater_than("external_travel_time", ett, 0, True)
        self._external_travel_time = float(ett)

    def to_xml_element(self):
        """Create an XML element for DNP drift settings.

        The velocity field data is written separately as a <Field> element.
        This element only references it by ID.

        Returns
        -------
        element : xml.etree._Element
            XML element representing DNP drift settings.

        """
        element = ET.Element("dnp_drift")

        vel_elem = ET.SubElement(element, "velocity_field")
        vel_elem.text = str(self.velocity_field.id)

        bm_elem = ET.SubElement(element, "boundary_map")
        for btype, surface_ids in self.boundary_map.items():
            bt_elem = ET.SubElement(bm_elem, btype)
            bt_elem.text = ' '.join(str(s) for s in surface_ids)

        mapping_elem = ET.SubElement(element, "mapping")
        mapping_elem.text = self.mapping

        integrator_elem = ET.SubElement(element, "integrator")
        integrator_elem.text = self.integrator

        dt_elem = ET.SubElement(element, "integrator_dt")
        dt_elem.text = str(self.integrator_dt)

        recycling_elem = ET.SubElement(element, "recycling")
        recycling_elem.text = str(self.recycling).lower()

        ett_elem = ET.SubElement(element, "external_travel_time")
        ett_elem.text = str(self.external_travel_time)

        return element

    @classmethod
    def from_xml_element(cls, elem, fields):
        """Construct a DNPDrift instance from an XML element.

        Parameters
        ----------
        elem : xml.etree._Element
            <dnp_drift> XML element to parse.
        fields : dict
            Dictionary mapping field IDs (int) to Field instances.

        Returns
        -------
        DNPDrift
            Reconstructed DNP drift settings.

        Raises
        ------
        ValueError
            If the referenced velocity field ID is not found or is
            not a VelocityField. Also, if the boundary_map is not found.

        """
        vel_elem = elem.find("velocity_field")
        if vel_elem is None:
            raise ValueError("Missing required <velocity_field> in <dnp_drift>.")
        vel_id = int(vel_elem.text)
        if vel_id not in fields:
            raise ValueError(
                f"DNP drift references velocity field id={vel_id}, "
                f"but no <Field> element with that id was found."
            )
        vel_field = fields[vel_id]
        if not isinstance(vel_field, VelocityField):
            raise ValueError(
                f"Field id={vel_id} referenced by <dnp_drift> is of "
                f"type '{vel_field._field_type}', expected 'velocity'."
            )

        kwargs = {"velocity_field": vel_field}

        bm = elem.find("boundary_map")
        if bm is None:
            raise ValueError("Missing required <boundary_map> in <dnp_drift>.")
        boundary_map = {}
        for child in bm:
            ids = [int(x) for x in child.text.split()]
            boundary_map[child.tag] = ids
        kwargs["boundary_map"] = boundary_map

        mapping = elem.find("mapping")
        if mapping is not None:
            kwargs["mapping"] = mapping.text

        integrator = elem.find("integrator")
        if integrator is not None:
            kwargs["integrator"] = integrator.text

        dt = elem.find("integrator_dt")
        if dt is not None:
            kwargs["integrator_dt"] = float(dt.text)

        recycling = elem.find("recycling")
        if recycling is not None:
            kwargs["recycling"] = recycling.text.lower() == "true"

        ett = elem.find("external_travel_time")
        if ett is not None:
            kwargs["external_travel_time"] = float(ett.text)

        return cls(**kwargs)

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        return (
            self.velocity_field == other.velocity_field
            and self.boundary_map == other.boundary_map
            and self.mapping == other.mapping
            and self.integrator == other.integrator
            and self.integrator_dt == other.integrator_dt
            and self.recycling == other.recycling
            and self.external_travel_time == other.external_travel_time
        )

    def __repr__(self):
        return (
            f"DNPDrift(velocity_field_id={self.velocity_field.id}, "
            f"boundary_map={self.boundary_map}, "
            f"mapping='{self.mapping}', "
            f"integrator='{self.integrator}', "
            f"integrator_dt={self.integrator_dt}, "
            f"recycling={self.recycling}, "
            f"external_travel_time={self.external_travel_time})"
        )
