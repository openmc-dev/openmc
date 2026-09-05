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
    physical_group_map : dict
        Mapping from mesh face IDs to physical group IDs. Must contain
        keys "face_ids" and "physical_groups", each mapping to a list
        of int of equal length.
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
    physical_group_map : dict
        Mapping from mesh face IDs to physical group IDs.
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
        physical_group_map,
        integrator="RK4",
        integrator_dt=0.1,
        recycling=True,
        external_travel_time=1.0,
    ):
        self.velocity_field = velocity_field
        self.boundary_map = boundary_map
        self.physical_group_map = physical_group_map
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
    def physical_group_map(self) -> dict[str, list[int]]:
        return self._physical_group_map

    @physical_group_map.setter
    def physical_group_map(self, physical_group_map: dict[str, list[int]]):
        cv.check_type("physical_group_map", physical_group_map, dict)

        if "face_ids" not in physical_group_map:
            raise ValueError("physical_group_map must contain 'face_ids'.")
        if "physical_groups" not in physical_group_map:
            raise ValueError("physical_group_map must contain 'physical_groups'.")

        face_ids = list(physical_group_map["face_ids"])
        physical_groups = list(physical_group_map["physical_groups"])

        for fid in face_ids:
            cv.check_type("face if in physical_group_map", fid, int)
        for pg in physical_groups:
            cv.check_type("physical group in physical_group_map", pg, int)

        if len(face_ids) != len(physical_groups):
            raise ValueError(
                f"physical_group_map: face_ids has {len(face_ids)} entries but "
                f"physical_groups has {len(physical_groups)} entries. "
                f"They must be identical."
            )

        self._physical_group_map = {
            "face_ids": face_ids,
            "physical_groups": physical_groups,
        }

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

        pgm_elem = ET.SubElement(element, "physical_group_map")
        face_ids_elem = ET.SubElement(pgm_elem, "face_ids")
        face_ids_elem.text = " ".join(
            str(i) for i in self.physical_group_map["face_ids"]
        )
        pg_elem = ET.SubElement(pgm_elem, "physical_groups")
        pg_elem.text = " ".join(
            str(i) for i in self.physical_group_map["physical_groups"]
        )

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

        # Physical group map (required)
        pgm_elem = elem.find("physical_group_map")
        if pgm_elem is None:
            raise ValueError("dnp_drift: missing <physical_group_map> element.")
        face_ids_elem = pgm_elem.find("face_ids")
        if face_ids_elem is None or face_ids_elem.text is None:
            raise ValueError(
                "dnp_drift: missing or empty <face_ids> in the <physical_group_map> element."
            )
        pg_elem = pgm_elem.find("physical_groups")
        if pg_elem is None or pg_elem.text is None:
            raise ValueError(
                "dnp_drift: missing or empty <physical_groups> in the <physical_group_map> element."
            )
        kwargs["physical_group_map"] = {
            "face_ids": [int(i) for i in face_ids_elem.text.split()],
            "physical_groups": [int(i) for i in pg_elem.text.split()]
        }

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
            and self.physical_group_map == other.physical_group_map
            and self.integrator == other.integrator
            and self.integrator_dt == other.integrator_dt
            and self.recycling == other.recycling
            and self.external_travel_time == other.external_travel_time
        )

    def __repr__(self):
        return (
            f"DNPDrift(velocity_field_id={self.velocity_field.id}, "
            f"boundary_map={self.boundary_map}, "
            f"physical_group_map=({len(self._physical_group_map['face_ids'])} faces), "
            f"integrator='{self.integrator}', "
            f"integrator_dt={self.integrator_dt}, "
            f"recycling={self.recycling}, "
            f"external_travel_time={self.external_travel_time})"
        )
