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

Verifications scripts are made available in **_verification**. They are designed
to compare the results of the regression tests to a strictly equivalent model
where the tracks of particles are expected to be identical. These scripts
should be run using similar settings than what is used with the regression tests.

.. warning:: The nested cubes verification scripts for dagmc are currently missing
    because there are errors with the geometry. These errors are currently under
    investigation and the scripts will be made available once the problem is fixed.
