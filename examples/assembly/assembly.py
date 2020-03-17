"""
This script builds a single PWR assembly and is a slightly more advanced
demonstration of model building using Python. The creation of two universes for
fuel pins and guide tube pins has been separated into functions, and then the
overall model is built by an `assembly` function. This script also demonstrates
the use of the `Model` class, which provides some extra convenience over using
`Geometry`, `Materials`, and `Settings` classes directly. Finally, the script
takes two command-line flags that indicate whether to build and/or run the
model.

"""

import argparse
from math import log10

import numpy as np
import openmc

# Define surfaces
fuel_or = openmc.ZCylinder(r=0.39218, name='Fuel OR')
clad_or = openmc.ZCylinder(r=0.45720, name='Clad OR')

# Define materials
fuel = openmc.Material(name='Fuel')
fuel.set_density('g/cm3', 10.29769)
fuel.add_nuclide('U234', 4.4843e-6)
fuel.add_nuclide('U235', 5.5815e-4)
fuel.add_nuclide('U238', 2.2408e-2)
fuel.add_nuclide('O16', 4.5829e-2)

clad = openmc.Material(name='Cladding')
clad.set_density('g/cm3', 6.55)
clad.add_nuclide('Zr90', 2.1827e-2)
clad.add_nuclide('Zr91', 4.7600e-3)
clad.add_nuclide('Zr92', 7.2758e-3)
clad.add_nuclide('Zr94', 7.3734e-3)
clad.add_nuclide('Zr96', 1.1879e-3)

hot_water = openmc.Material(name='Hot borated water')
hot_water.set_density('g/cm3', 0.740582)
hot_water.add_nuclide('H1', 4.9457e-2)
hot_water.add_nuclide('O16', 2.4672e-2)
hot_water.add_nuclide('B10', 8.0042e-6)
hot_water.add_nuclide('B11', 3.2218e-5)
hot_water.add_s_alpha_beta('c_H_in_H2O')


def fuel_pin():
    """Returns a fuel pin universe."""

    fuel_cell = openmc.Cell(fill=fuel, region=-fuel_or)
    clad_cell = openmc.Cell(fill=clad, region=+fuel_or & -clad_or)
    hot_water_cell = openmc.Cell(fill=hot_water, region=+clad_or)

    univ = openmc.Universe(name='Fuel Pin')
    univ.add_cells([fuel_cell, clad_cell, hot_water_cell])
    return univ


def guide_tube_pin():
    """Returns a control rod guide tube universe."""

    gt_inner_cell = openmc.Cell(fill=hot_water, region=-fuel_or)
    gt_clad_cell = openmc.Cell(fill=clad, region=+fuel_or & -clad_or)
    gt_outer_cell = openmc.Cell(fill=hot_water, region=+clad_or)

    univ = openmc.Universe(name='Guide Tube')
    univ.add_cells([gt_inner_cell, gt_clad_cell, gt_outer_cell])
    return univ


def assembly_model():
    """Returns a single PWR fuel assembly."""

    model = openmc.model.Model()

    # Create fuel assembly Lattice
    pitch = 21.42
    assembly = openmc.RectLattice(name='Fuel Assembly')
    assembly.pitch = (pitch/17, pitch/17)
    assembly.lower_left = (-pitch/2, -pitch/2)

    # Create array indices for guide tube locations in lattice
    gt_pos = np.array([
                  [2, 5],   [2, 8],   [2, 11],
             [3, 3],                        [3, 13],
        [5, 2],   [5, 5],   [5, 8],   [5, 11],    [5, 14],
        [8, 2],   [8, 5],   [8, 8],   [8, 11],    [8, 14],
        [11, 2],  [11, 5],  [11, 8],  [11, 11],   [11, 14],
             [13, 3],                       [13, 13],
                  [14, 5],  [14, 8],  [14, 11]
    ])

    # Create 17x17 array of universes. First we create a 17x17 array all filled
    # with the fuel pin universe. Then, we replace the guide tube positions with
    # the guide tube pin universe (note the use of numpy fancy indexing to
    # achieve this).
    assembly.universes = np.full((17, 17), fuel_pin())
    assembly.universes[gt_pos[:, 0], gt_pos[:, 1]] = guide_tube_pin()

    # Create outer boundary of the geometry to surround the lattice
    outer_boundary = openmc.model.rectangular_prism(
        pitch, pitch, boundary_type='reflective')

    # Create a cell filled with the lattice
    main_cell = openmc.Cell(fill=assembly, region=outer_boundary)

    # Finally, create geometry by providing a list of cells that fill the root
    # universe
    model.geometry = openmc.Geometry([main_cell])

    model.settings.batches = 150
    model.settings.inactive = 50
    model.settings.particles = 1000
    model.settings.source = openmc.Source(space=openmc.stats.Box(
        (-pitch/2, -pitch/2, -1),
        (pitch/2, pitch/2, 1),
        only_fissionable=True
    ))

    # NOTE: We never actually created a Materials object. When you export/run
    # using the Model object, if no materials were assigned it will look through
    # the Geometry object and automatically export any materials that are
    # necessary to build the model.
    return model


if __name__ == '__main__':
    # Set up command-line arguments for generating/running the model
    parser = argparse.ArgumentParser()
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()

    if args.generate or args.run:
        model = assembly_model()
        if args.generate:
            model.export_to_xml()
        if args.run:
            model.run()
