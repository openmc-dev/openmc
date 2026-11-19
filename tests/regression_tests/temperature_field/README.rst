Temperature field regression tests
==================================

These regression tests use four models selected to cover distinct interactions
without testing the full product of geometry representations and boundary
conditions:

* A CSG cube with reflective boundaries tests coincident geometry and
  temperature-field boundaries.
* Nested CSG cubes with periodic boundaries test entry into and exit from the
  temperature field, cell-temperature fallback, and relocation after a
  periodic translation.
* Nested cubes represented by a lattice with vacuum boundaries test lattice
  crossings and particle leakage.
* Nested cubes represented by DAGMC with reflective boundaries test DAGMC ray
  history and reflection.

The temperature field uses a 2x2x2 regular mesh where each cell has its own
temperature value (built by starting with 294K and adding increments of 100K
to assign the next cell's temperature). When modeled, the temperature of the
surrounding cube is set to 250K.

For the DAGMC geometry, the boundary conditions and temperature values are
directly stored in the DAGMC file created with a modified version of
StellarMesh.

Independent comparisons against explicit cell-temperature models are integration
tests in ``tests/unit_tests/test_temperature_field_verification.py``. They
compare the two configurations directly, without golden files, so regenerating
these regressions cannot change the expected equivalence. Run those integration
tests when updating regression references.
