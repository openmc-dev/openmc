Temperature-field verification assets
=====================================

``nested_cubes_multicell_reflecting.h5m`` is the explicit cell-temperature
reference geometry from Joffrey Dorville's verification repository, commit
``80937e2c740071fcf5fb428c1fa8551624ed9801``, at
``openmc_models/dagmc/nested_cubes/reflective_bc``:

https://github.com/JoffreyDorville/temperature_field_verification

The shared builders in ``tests/temperature_field.py`` adapt that repository's
CSG and lattice reference models. The integration tests compare eigenvalue,
tally sums, and tally sums of squares directly, allowing only floating-point
rounding (relative tolerance 1e-12). They use no stored result baselines and
run with the executable and history/event mode selected by pytest.

The asset and adapted reference models are distributed under the accompanying
MIT license. The mesh is checked in so tests need no CAD-generation packages.
