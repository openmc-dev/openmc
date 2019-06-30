"""
Helper class for building a pin from concentric cylinders
"""


from math import sqrt
from numbers import Real
from operator import attrgetter

import openmc
import openmc.checkvalue as cv
from . import subdivide

__all__ = ["Pin"]


class Pin(openmc.Universe):
    """Special universe used to model pins

    Parameters
    ----------
    surfaces: iterable of :class:`openmc.Cylinder`
        Cylinders used to define boundaries
        between materials. All cylinders must be
        concentric and of the same orientation, e.g.
        all :class:`openmc.ZCylinders`
    materials: iterable of :class:`openmc.Material`
        Materials to go between ``surfaces``. There must be one
        more material than surfaces, corresponding to the material
        that spans all space outside the final ring.
    universe_id: None or int
        Unique identifier for this universe
    name: str
        Name for this universe

    See Also
    --------

    :meth:`openmc.Pin.from_radii` - Convinience function to build
    from radii of cylinders
    """

    def __init__(self, surfaces, materials, universe_id=None, name=""):
        cv.check_iterable_type("materials", materials, openmc.Material)
        cv.check_length("surfaces", surfaces, len(materials) - 1)

        # Ensure that all surfaces are same type of cylinder
        self._check_surfaces(surfaces)
        regions = subdivide(surfaces)

        cells = [
            openmc.Cell(fill=m, region=r) for m, r in zip(materials, regions)
        ]

        super().__init__(universe_id=universe_id, name=name, cells=cells)
        self._list_cells = cells  # need ordered by radial position

    def _check_surfaces(self, surfaces):
        cv.check_type(
            "surface 0",
            surfaces[0],
            (openmc.ZCylinder, openmc.YCylinder, openmc.XCylinder),
        )
        if isinstance(surfaces[0], openmc.ZCylinder):
            center_getter = attrgetter("x0", "y0")
        elif isinstance(surfaces[0], openmc.YClylinder):
            center_getter = attrgetter("x0", "y0")
        elif isinstance(surfaces[0], openmc.XClylinder):
            center_getter = attrgetter("z0", "y0")
        else:
            raise TypeError(
                "Not configured to interpret {} surfaces".format(
                    surfaces[0].__class__.__name__
                )
            )
        cv.check_iterable_type("surfaces", surfaces[1:], type(surfaces[0]))
        # Check for concentric-ness and increasing radii
        centers = set()
        radii = []
        rad = 0
        for ix, surf in enumerate(surfaces):
            cur_rad = surf.r
            if cur_rad <= rad:
                raise ValueError(
                    "Surfaces do not appear to be increasing in radius. "
                    "Surface {} at index {} has radius {:7.3E} compared to "
                    "previous radius of {:7.5E}".format(
                        surf.id, ix, cur_rad, rad
                    )
                )
            rad = cur_rad
            radii.append(cur_rad)
            centers.add(center_getter(surf))

        if len(centers) > 1:
            raise ValueError(
                "Surfaces do not appear to be concentric. The following "
                "centers were found: {}".format(centers)
            )
        self._radii = radii
        self._surfaces = surfaces
        self._surf_type = type(surfaces[0])

    @classmethod
    def from_radii(
        cls,
        radii,
        materials,
        universe_id=None,
        name="",
        orientation="z",
        center=(0.0, 0.0),
    ):
        """Construct using radii of concentric cylinders and materials

        Parameters
        ----------
        radii: iterable of float
            Radii of the intended cylinders. Must be all positive
            values and increasing
        materials: iterable of :class:`openmc.Material`
            Materials used to fill cylinders created by ``radii``,
            starting from inside to the outside of the pin.
            There must be one extra material corresponding to all
            area outside the last ring
        universe_id: None or int
            Unique identifier to give this pin
        name: str
            Name of the created universe
        orientation: {"x", "y", "z"}
            Axis along which to orient the pin. Default is ``"z"``
        center: iterable of float
            Center of the pin in the plane perpendicular
        """
        cv.check_iterable_type("materials", materials, openmc.Material)
        cv.check_length("radii", radii, len(materials) - 1)
        for ix, rad in enumerate(radii):
            if rad < 0:
                raise ValueError(
                    "Radius {:7.3E} at index {} is non-positive".format(
                        rad, ix
                    )
                )
            if ix and rad <= radii[ix - 1]:
                raise ValueError(
                    "Radii must be increasing values. Radius {:7.3E} at index "
                    "{} is not greater than previous value of {:7.3E}".format(
                        rad, ix, radii[ix - 1]
                    )
                )

        if orientation == "z":
            surfCls = openmc.ZCylinder
            basis = ["x0", "y0"]
        elif orientation == "y":
            surfCls = openmc.YCylinder
            basis = ["x0", "z0"]
        elif orientation == "z":
            surfCls = openmc.XCylinder
            basis = ["y0", "z0"]
        else:
            raise ValueError(
                "Orientation of {} not understood".format(orientation)
            )

        centerKwargs = dict(zip(basis, center))

        surfaces = [surfCls(r=rad, **centerKwargs) for rad in radii]

        return cls(surfaces, materials, universe_id=universe_id, name=name)

    def subdivide_ring(self, ring_index, n_divs):
        """Divide one ring of the pin into equal-area rings

        Each new ring will be added to the model, and filled with
        a unique material copied from the original.

        Parameters
        ring_index: int
            Index of the ring to be divided where 0 is the innermost
            ring. Will not divide the outermost region as there is no
            upper bound
        n_divs: int
            Number of equal area divisions to make in this ring
        """
        # Don't allow subdivision of outer, infinite region
        cv.check_less_than("ring_index", ring_index, len(self._radii))
        cv.check_type("n_divs", n_divs, Real)
        cv.check_greater_than("n_divs", n_divs, 1)

        if ring_index < 0:
            ring_index = len(self._list_cells) + ring_index

        # Get all the information we need to replicate this
        # region with unique insides
        orig_cell = self._list_cells[ring_index]

        lower_rad = self._radii[ring_index - 1] if ring_index else 0.0
        area_term = (self._radii[ring_index] ** 2 - lower_rad ** 2) / n_divs

        new_radii = []
        new_surfaces = []
        new_cells = []

        # Adding N - 1 new regions
        # N - 2 surfaces are made
        # Original cell is not removed, but not occupies last ring

        for i in range(n_divs - 1):
            r = sqrt(area_term + lower_rad ** 2)
            lower_rad = r
            new_radii.append(r)
            surf = self._surf_type(r=r)
            new_surfaces.append(surf)
            if i == 0:
                if ring_index:
                    region = (
                        -surf & +self._surfaces[ring_index - 1]
                    )
                else:
                    region = -surf
            else:
                region = -surf & +new_surfaces[-2]
            new_cells.append(
                openmc.Cell(region=region, fill=orig_cell.fill.clone())
            )

        orig_cell.region = -self._surfaces[ring_index] & +surf

        self.add_cells(new_cells)

        self._list_cells = (
            self._list_cells[:ring_index]
            + new_cells
            + self._list_cells[ring_index:]
        )
        self._radii = (
            self._radii[:ring_index] + new_radii + self._radii[ring_index:]
        )
        self._surfaces = (
            self._surfaces[:ring_index]
            + new_surfaces
            + self._surfaces[ring_index:]
        )

    @property
    def radii(self):
        """Return a tuple of the radii in this :class:`Pin`"""
        return tuple(self._radii)
