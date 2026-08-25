from collections.abc import Iterable
from math import sqrt
from itertools import product
import numpy as np
from operator import attrgetter
from warnings import warn

from openmc import Cylinder, Universe, Cell
from .surface_composite import RectangularPrism, HexagonalPrism
from ..checkvalue import (check_type, check_value, check_length,
                          check_less_than, check_iterable_type)
import openmc.data


ZERO_CELSIUS_TO_KELVIN = 273.15
ZERO_FAHRENHEIT_TO_KELVIN = 459.67
PSI_TO_MPA = 0.006895


def borated_water(boron_ppm, temperature=293., pressure=0.1013, temp_unit='K',
                  press_unit='MPa', density=None, **kwargs):
    """Return a Material with the composition of boron dissolved in water.

    The water density can be determined from a temperature and pressure, or it
    can be set directly.

    The concentration of boron has no effect on the stoichiometric ratio of H
    and O---they are fixed at 2-1.

    Parameters
    ----------
    boron_ppm : float
        The weight fraction in parts-per-million of elemental boron in the
        water.
    temperature : float
        Temperature in [K] used to compute water density.
    pressure : float
        Pressure in [MPa] used to compute water density.
    temp_unit : {'K', 'C', 'F'}
        The units used for the `temperature` argument.
    press_unit : {'MPa', 'psi'}
        The units used for the `pressure` argument.
    density : float
        Water density in [g / cm^3].  If specified, this value overrides 
        the value that is computed from the temperature and pressure arguments.
    **kwargs
        All keyword arguments are passed to the created Material object.

    Returns
    -------
    openmc.Material

    """
    # Perform any necessary unit conversions.
    check_value('temperature unit', temp_unit, ('K', 'C', 'F'))
    if temp_unit == 'K':
        T = temperature
    elif temp_unit == 'C':
        T = temperature + ZERO_CELSIUS_TO_KELVIN
    elif temp_unit == 'F':
        T = (temperature + ZERO_FAHRENHEIT_TO_KELVIN) * (5/9)
    check_value('pressure unit', press_unit, ('MPa', 'psi'))
    if press_unit == 'MPa':
        P = pressure
    elif press_unit == 'psi':
        P = pressure * PSI_TO_MPA

    # Set the density of water, either from an explicitly given density or from
    # temperature and pressure.
    if density is not None:
        water_density = density
    else:
        water_density = openmc.data.water_density(T, P)

    # Compute the density of the solution.
    solution_density = water_density / (1 - boron_ppm * 1e-6)

    # Compute the molar mass of pure water.
    hydrogen = openmc.Element('H')
    oxygen = openmc.Element('O')
    M_H2O = 0.0
    for iso_name, frac, junk in hydrogen.expand(2.0, 'ao'):
        M_H2O += frac * openmc.data.atomic_mass(iso_name)
    for iso_name, frac, junk in oxygen.expand(1.0, 'ao'):
        M_H2O += frac * openmc.data.atomic_mass(iso_name)

    # Compute the molar mass of boron.
    boron = openmc.Element('B')
    M_B = 0.0
    for iso_name, frac, junk in boron.expand(1.0, 'ao'):
        M_B += frac * openmc.data.atomic_mass(iso_name)

    # Compute the number fractions of each element.
    frac_H2O = (1 - boron_ppm * 1e-6) / M_H2O
    frac_H = 2 * frac_H2O
    frac_O = frac_H2O
    frac_B = boron_ppm * 1e-6 / M_B

    # Build the material.
    out = openmc.Material(temperature=T, **kwargs)
    out.add_element('H', frac_H, 'ao')
    out.add_element('O', frac_O, 'ao')
    out.add_element('B', frac_B, 'ao')
    out.set_density('g/cc', solution_density)
    out.add_s_alpha_beta('c_H_in_H2O')
    return out




def rectangular_prism(width, height, axis='z', origin=(0., 0.),
                      boundary_type='transmission', corner_radius=0.):
    warn("The rectangular_prism(...) function has been replaced by the "
         "RectangularPrism(...) class. Future versions of OpenMC will not "
         "accept rectangular_prism.", FutureWarning)
    return -RectangularPrism(
        width=width, height=height, axis=axis, origin=origin,
        boundary_type=boundary_type, corner_radius=corner_radius)


def hexagonal_prism(edge_length=1., orientation='y', origin=(0., 0.),
                    boundary_type='transmission', corner_radius=0.):
    warn("The hexagonal_prism(...) function has been replaced by the "
         "HexagonalPrism(...) class. Future versions of OpenMC will not "
         "accept hexagonal_prism.", FutureWarning)
    return -HexagonalPrism(
        edge_length=edge_length, orientation=orientation, origin=origin,
        boundary_type=boundary_type, corner_radius=corner_radius)


def get_hexagonal_prism(*args, **kwargs):
    warn("get_hexagonal_prism(...) has been renamed hexagonal_prism(...). "
         "Future versions of OpenMC will not accept get_hexagonal_prism.",
         FutureWarning)
    return hexagonal_prism(*args, **kwargs)


cylinder_from_points = Cylinder.from_points


def subdivide(surfaces):
    """Create regions separated by a series of surfaces.

    This function allows regions to be constructed from a set of a surfaces that
    are "in order". For example, if you had four instances of
    :class:`openmc.ZPlane` at z=-10, z=-5, z=5, and z=10, this function would
    return a list of regions corresponding to z < -10, -10 < z < -5, -5 < z < 5,
    5 < z < 10, and 10 < z. That is, for n surfaces, n+1 regions are returned.

    Parameters
    ----------
    surfaces : sequence of openmc.Surface
        Surfaces separating regions

    Returns
    -------
    list of openmc.Region
        Regions formed by the given surfaces

    """
    regions = [-surfaces[0]]
    for s0, s1 in zip(surfaces[:-1], surfaces[1:]):
        regions.append(+s0 & -s1)
    regions.append(+surfaces[-1])
    return regions


def pin(surfaces, items, subdivisions=None, divide_vols=True,
        **kwargs):
    """Convenience function for building a fuel pin

    Parameters
    ----------
    surfaces : iterable of :class:`openmc.Cylinder`
        Cylinders used to define boundaries
        between items. All cylinders must be
        concentric and of the same orientation, e.g.
        all :class:`openmc.ZCylinder`
    items : iterable
        Objects to go between ``surfaces``. These can be anything
        that can fill a :class:`openmc.Cell`, including
        :class:`openmc.Material`, or other :class:`openmc.Universe`
        objects. There must be one more item than surfaces,
        which will span all space outside the final ring.
    subdivisions : None or dict of int to int
        Dictionary describing which rings to subdivide and how
        many times. Keys are indexes of the annular rings
        to be divided. Will construct equal area rings
    divide_vols : bool
        If this evaluates to ``True``, then volumes of subdivided
        :class:`openmc.Material` instances will also be divided by the
        number of divisions.  Otherwise the volume of the
        original material will not be modified before subdivision
    kwargs:
        Additional key-word arguments to be passed to
        :class:`openmc.Universe`, like ``name="Fuel pin"``

    Returns
    -------
    :class:`openmc.Universe`
        Universe of concentric cylinders filled with the desired
        items
    """
    if "cells" in kwargs:
        raise ValueError(
            "Cells will be set by this function, not from input arguments.")
    check_type("items",  items, Iterable)
    check_length("surfaces", surfaces, len(items) - 1, len(items) - 1)
    # Check that all surfaces are of similar orientation
    check_type("surface", surfaces[0], Cylinder)
    surf_type = type(surfaces[0])
    check_iterable_type("surfaces", surfaces[1:], surf_type)

    # Check for increasing radii and equal centers
    if surf_type is openmc.ZCylinder:
        center_getter = attrgetter("x0", "y0")
    elif surf_type is openmc.YCylinder:
        center_getter = attrgetter("x0", "z0")
    elif surf_type is openmc.XCylinder:
        center_getter = attrgetter("z0", "y0")
    else:
        raise TypeError(
            f"Not configured to interpret {surf_type.__name__} surfaces")

    centers = set()
    prev_rad = 0
    for ix, surf in enumerate(surfaces):
        cur_rad = surf.r
        if cur_rad <= prev_rad:
            raise ValueError(
                "Surfaces do not appear to be increasing in radius. "
                "Surface {} at index {} has radius {:7.3e} compared to "
                "previous radius of {:7.5e}".format(
                    surf.id, ix, cur_rad, prev_rad))
        prev_rad = cur_rad
        centers.add(center_getter(surf))

    if len(centers) > 1:
        raise ValueError(
            "Surfaces do not appear to be concentric. The following "
            "centers were found: {}".format(centers))

    if subdivisions is not None:
        check_length("subdivisions", subdivisions, 1, len(surfaces))
        orig_indexes = list(subdivisions.keys())
        check_iterable_type("ring indexes", orig_indexes, int)
        check_iterable_type(
            "number of divisions", list(subdivisions.values()), int)
        for ix in orig_indexes:
            if ix < 0:
                subdivisions[len(surfaces) + ix] = subdivisions.pop(ix)
        # Dissallow subdivision on outer most, infinite region
        check_less_than(
            "outer ring", max(subdivisions), len(surfaces), equality=True)

        # ensure ability to concatenate
        if not isinstance(items, list):
            items = list(items)
        if not isinstance(surfaces, list):
            surfaces = list(surfaces)

        # generate equal area divisions
        # Adding N - 1 new regions
        # N - 2 surfaces are made
        # Original cell is not removed, but not occupies last ring
        for ring_index in reversed(sorted(subdivisions.keys())):
            nr = subdivisions[ring_index]
            new_surfs = []

            lower_rad = 0.0 if ring_index == 0 else surfaces[ring_index - 1].r

            upper_rad = surfaces[ring_index].r

            area_term = (upper_rad ** 2 - lower_rad ** 2) / nr

            for new_index in range(nr - 1):
                lower_rad = sqrt(area_term + lower_rad ** 2)
                new_surfs.append(surf_type(r=lower_rad))

            surfaces = (
                    surfaces[:ring_index] + new_surfs + surfaces[ring_index:])

            filler = items[ring_index]
            if (divide_vols and hasattr(filler, "volume")
                    and filler.volume is not None):
                filler.volume /= nr

            items[ring_index:ring_index] = [
                filler.clone() for _i in range(nr - 1)]

    # Build the universe
    regions = subdivide(surfaces)
    cells = [Cell(fill=f, region=r) for r, f in zip(regions, items)]
    return Universe(cells=cells, **kwargs)

# Pointset generators for discrete-ordinate angular meshes

def levelsymmetric_sn(N: int, mu1_sq: float = None):
    """
    Generates the direction vectors of an order-N level-symmetric 
    quadrature set. Does not provide the corresponding weights.

    Parameters
    ----------
    N : int
        Even quadrature order (e.g. 2, 4, 6, ..., 20).
    mu1_sq : float
        Square of the first direction cosine mu_1, used to generate other 
        levels. Must lie in the open interval (0, 1/3). For the special case 
        N == 2, mu1_sq must be equal to 1/3.

    Returns
    -------
    np.ndarray
        Array of shape (N*(N+2), 3) giving the [x, y, z] unit-vector
        directions of every point in the quadrature set.
    """
    if N < 2 or N % 2 != 0:
        raise ValueError(f"N must be a positive even integer, got {N}")
    M = N // 2

    # lookup mu1_sq if not provided
    if mu1_sq is not None:
        if N == 2:
            if not np.isclose(mu1_sq, 1 / 3):
                raise ValueError("For N=2, mu1_sq must equal 1/3 (got {mu1_sq}).")
            mu = np.array([np.sqrt(1 / 3)])
        elif N > 20:
            raise ValueError(
                f"Level-symmetric quadrature generation only supported for 2<=N<=20; got {N}")
        else:
            if not (0 < mu1_sq < 1 / 3):
                raise ValueError(
                    f"mu1_sq={mu1_sq} is out of the valid range (0, 1/3) for N={N}."
                )
            
            Delta = (1 - 3 * mu1_sq) / (M - 1)
            mu_sq = np.array([mu1_sq + i * Delta for i in range(M)])
            mu = np.sqrt(mu_sq)

    else:
        match N:
            case 2:
                pass
            case 4:
                mu1_sq = 0.1225148226554413
            case 6:
                mu1_sq = 0.0710944373419735
            case 8:
                mu1_sq = 0.0476190476190470
            case 10:
                mu1_sq = 0.0358425646593916
            case 12:
                mu1_sq = 0.0279600712640057
            case 14:
                mu1_sq = 0.0230997020840970
            case 16:
                mu1_sq = 0.0193090131285642
            case 18:
                mu1_sq = 0.0167300008552435
            case 20:
                mu1_sq = 0.0145451663522475
            case _:
                raise ValueError(
                    f"Level-symmetric quadrature generation only supported for 2<=N<=20; got {N}")

        if M == 1:
            mu = np.array([np.sqrt(1 / 3)])
        else:
            Delta = (1 - 3 * mu1_sq) / (M - 1)
            mu_sq = np.array([mu1_sq + i * Delta for i in range(M)])
            mu = np.sqrt(mu_sq)

    # generate the angles for 1 octant of unit sphere
    octant_points = []
    for l in range(0, M):
        for m in range(0, M):
            for n in range(0, M):
                if l + m + n == (M + 1):
                    octant_points.append((mu[l], mu[m], mu[n]))
    octant_points = np.array(octant_points)

    # reflect into other octants
    signs = list(product([1, -1], repeat=3))
    all_points = np.array(
        [pt * np.array(s) for pt in octant_points for s in signs]
    )

    return all_points

def tcl_sn(N: int):
    """
    Generate the direction vectors of an order-N triangular 
    Chebyshev-Legendre (TCL) quadrature set. Supports even N >= 4.

    Parameters
    ----------
    N : int
        Even quadrature order (e.g. 4, 6, 8, ...).

    Returns
    -------
    np.ndarray
        Array of shape (N*(N+2), 3) giving the [x, y, z] unit-vector
        directions of every point in the quadrature set.
    """
    if N < 4 or N % 2 != 0:
        raise ValueError(f"N must be an even integer >= 4, got {N}")

    M = N // 2

    # get polar levels
    nodes, _ = np.polynomial.legendre.leggauss(N)
    mu = np.sort(nodes[nodes > 0])

    # generate the angles for 1 octant using Chebyshev quadrature for 
    # azimuthal angles
    octant_points = []
    for i in range(M):
        count = M - i # rings get smaller towards pole
        m = mu[i]
        sin_theta = np.sqrt(1 - m**2)
        for k in range(1, count + 1):
            phi = (2 * k - 1) * (np.pi / 2) / (2 * count)
            x = sin_theta * np.cos(phi)
            y = sin_theta * np.sin(phi)
            z = m
            octant_points.append((x, y, z))
    octant_points = np.array(octant_points)

    # reflect into other octants
    signs = list(product([1, -1], repeat=3))
    all_points = np.array(
        [pt * np.array(s) for pt in octant_points for s in signs]
    )

    assert np.allclose(np.linalg.norm(all_points, axis=1), 1.0), (
        "Not all generated points lie on the unit sphere"
    )

    return all_points

def _subdivide_icosahedron_faces(vertices, faces, nu):
    """
    Given a list of the coordinates of the vertices of a unit icosahedron, 
    and a list linking sets of these vertices to individual faces of the unit 
    icosahedron, this function will subdivides each edge of the icosahedron 
    into nu equal segments, adding vertices on the edges and faces to produce 
    triangular subfaces of equal size.

    Parameters
    ----------
    vertices : numpy array of shape (n_verts, 3)
        Coordinates of each vertex of the unit icosahedron
    faces : numpy array of shape (n_faces, 3)
        List of vertex indices corresponding to each face of the base 
        icosahedron
    nu : int
        Subdivision frequency, integer > 1

    Returns
    -------
    subvertices : numpy array of shape (n_verts + n_faces*(nu+1)*(nu-1)/2, 3)
        List of vertices on subdivided icosahedron
    subfaces : numpy array of shape (n_faces*nu**2, 3)
        List of vertex indices corresponding to each face of the subdivided 
        icosahedron
    """
    edges = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [0, 2]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)

    n_faces = faces.shape[0]
    n_vertices = vertices.shape[0]
    n_edges = edges.shape[0]
    n_int_verts = (nu - 1) * (nu - 2) // 2

    n_subverts  = n_vertices + n_edges * (nu - 1) + n_faces * n_int_verts
    subvertices = np.empty((n_subverts, 3))
    subvertices[:n_vertices] = vertices

    # populate edge vertices:
    # position of the k-th vertex along edge AB is given by 
    # (1 − w_k)·a  +  w_k·b,   w_k = (k+1)/nu
    w  = np.arange(1, nu) / nu
    vA = vertices[edges[:, 0]]
    vB = vertices[edges[:, 1]]
    edge_verts = (1 - w)[:, None, None] * vA[None] + w[:, None, None] * vB[None]
    subvertices[n_vertices : n_vertices + n_edges * (nu - 1)] = (
        edge_verts.transpose(1, 0, 2).reshape(-1, 3)
    )

    f_A, f_B, f_C = faces[:, 0], faces[:, 1], faces[:, 2]

    edge_dict = {(int(a), int(b)):  i for i, (a, b) in enumerate(edges)}
    edge_dict.update({(int(b), int(a)): ~i for i, (a, b) in enumerate(edges)})

    def directed_edge_indices(u_arr, v_arr):
        # Return (n_faces, nu-1) global vertex indices along the u->v edge
        ei   = np.array([edge_dict[(int(u), int(v))] for u, v in zip(u_arr, v_arr)])
        base = n_vertices + np.where(ei >= 0, ei, ~ei) * (nu - 1)
        idx  = base[:, None] + np.arange(nu - 1)
        idx[ei < 0] = idx[ei < 0, ::-1]
        return idx

    AB = directed_edge_indices(f_A, f_B)
    AC = directed_edge_indices(f_A, f_C)
    BC = directed_edge_indices(f_B, f_C)

    # global indices of subvertices:
    # (0,0) = corner "A," (nu, 0) = corner B, (nu, nu) = corner C
    # and along the edges above
    local_idx = np.empty((n_faces, nu + 1, nu + 1), dtype=int)
    local_idx[:, 0,    0   ] = f_A
    local_idx[:, nu,   0   ] = f_B
    local_idx[:, nu,   nu  ] = f_C
    local_idx[:, 1:nu, 0   ] = AB 
    r_e = np.arange(1, nu)
    local_idx[:, r_e,  r_e ] = AC
    local_idx[:, nu,   1:nu] = BC

    # row, column indices of interior points
    r_int = np.array([r for r in range(2, nu) for _ in range(r - 1)])  # (n_int,)
    c_int = np.array([c for r in range(2, nu) for c in range(1, r)])
    T_base = n_vertices + n_edges * (nu - 1)

    if n_int_verts > 0:
        T_start = T_base + np.arange(n_faces) * n_int_verts
        local_idx[:, r_int, c_int] = (
            T_start[:, None] + np.arange(n_int_verts)[None, :]
        )

    # populate subfaces with vertex coordination info:
    # local vertex (r, c) with 0 ≤ c ≤ r ≤ nu has barycentric weights
    #   A'=(nu-r)/nu,  B'=(r-c)/nu,  C'=c/nu
    # interior vertices have 0 < col < row < nu
    tri_list = []
    for i in range(nu):
        for j in range(i):
            tri_list.append([(i,j), (i+1,j), (i+1,j+1)])
            tri_list.append([(i,j), (i+1,j+1), (i,j+1)])
        tri_list.append([(i,i), (i+1,i), (i+1,i+1)])
    tri_arr = np.array(tri_list)
    tr, tc  = tri_arr[:, :, 0], tri_arr[:, :, 1]
    subfaces = local_idx[:, tr, tc].reshape(n_faces * nu ** 2, 3)

    if n_int_verts > 0:
        alpha = (nu - r_int) / nu
        beta  = (r_int - c_int) / nu
        gamma = c_int / nu
        int_verts = (  alpha[None, :, None] * vertices[f_A][:, None, :] 
                     + beta [None, :, None] * vertices[f_B][:, None, :]
                     + gamma[None, :, None] * vertices[f_C][:, None, :])
        subvertices[T_base:] = int_verts.reshape(-1, 3)

    return subvertices, subfaces

def icosphere_sn(nu: int = 1, point_type: str = "centroids"):
    """
    Generates direction vectors from the vertices or centroids of a 
    unit spherical icosahedron with principal faces subdivided at the nu-th 
    frequency. 
    
    That is, beginning from a "parent" icosahedron inscribed in the unit 
    sphere, the edges are first divided into nu equal segments, and then the 
    endpoints of these segments are connected to subdivide each "parent" face
    into nu^2 triangular subfaces. Lastly, the vectors describing the 
    locations of either the vertices or centroids of these faces are 
    projected back onto the surface of the unit sphere.

    Parameters
    ----------
    nu : int
        Subdivision frequency
    point_type: str
        Keyword specifying whether to return the vertices ("vertices", 
        "vertex", or "vert") or centroids ("centroids", "centroid", or 
        "cent") of the subtriangles on the icosphere

    Returns
    -------
    pointset : numpy array of shape (12 + 10 * (nu+1) * (nu-1), 3) if
        specifying "vertices" or of shape (20 * nu**2, 3) if specifying 
        "centroids"

    """
    # check pointset type
    match point_type:
        case "vertices" | "vertex" | "vert" :
            return_type = "vert"
        case "centroids" | "centroid" | "cent":
            return_type = "cent"
        case _:
            return ValueError(f"Unknown pointset type specified (got {point_type})")

    # vertices of base icosahedron
    phi = (1 + np.sqrt(5)) / 2
    vertices = np.array([
        [0, 1, phi], [0, -1, phi], [1, phi, 0],
        [-1, phi, 0], [phi, 0, 1], [-phi, 0, 1],
        [0, -1, -phi], [0, 1, -phi], [-1, -phi, 0],
        [1, -phi, 0], [-phi, 0, -1], [phi, 0, -1]])
    vertices /= np.sqrt(1 + phi ** 2)
    
    # coordination of vertices to faces:
    # each entry in the array is a row vector, corresponding to one 
    # individual face of the icosahedron, giving the indices into the 
    # "vertices" array of its own particular vertices
    faces = np.array([
        [0, 5, 1], [0, 3, 5], [0, 2, 3], [0, 4, 2], [0, 1, 4],
        [1, 5, 8], [5, 3, 10], [3, 2, 7], [2, 4, 11], [4, 1, 9],
        [7, 11, 6], [11, 9, 6], [9, 8, 6], [8, 10, 6], [10, 7, 6],
        [2, 11, 7], [4, 9, 11], [1, 8, 9], [5, 10, 8], [3, 7, 10]])

    # subdividing
    if nu > 1:
        vertices, faces = _subdivide_icosahedron_faces(vertices, faces, nu)
        # project back to unit length
        vertices = vertices / np.sqrt(np.sum(vertices ** 2, axis=1, keepdims=True))
    
    if return_type == "vert":
        pointset = vertices
    else:
        # return centroids
        pointset = vertices[faces].mean(axis=1)
        pointset = pointset / np.sqrt(np.sum(pointset ** 2, axis=1, keepdims=True))

    return pointset