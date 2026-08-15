Temperature field regression tests
==================================

These regressions tests use two models:
- single_cube: a 5x5x5cm cube,
- nested_cubes: a 5x5x5cm cube surrounded by a 15x15x15cm cube.

Three types of boundary conditions are tested:
- vacuum,
- reflective,
- periodic (except with DAGMC geometries).

The temperature field uses a 2x2x2 regular mesh where each cell has its own
temperature value (built by starting with 294K and adding increments of 100K
to assign the next cell's temperature). When modeled, the temperature of the
surrounding cube is set to 250K.

For DAGMC geometries, the boundary conditions and temperature values are
directly stored in the DAGMC files created with a modified version of
stellarmesh.

The regression tests compare against reference input and result files generated
from equivalent temperature-field models. Use the standard regression-test
workflow to regenerate the reference files when the expected results change.
