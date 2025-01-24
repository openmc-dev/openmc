.. _variance_reduction:

==================
Variance Reduction
==================

Global variance reduction in OpenMC is accomplished by weight windowing
techniques. OpenMC is capable of generating weight windows using either the
MAGIC or FW-CADIS methods. Both techniques will produce a ``weight_windows.h5``
file that can be loaded and used later on. In this section, we break down the
steps required to both generate and then apply weight windows.

.. _ww_generator:

------------------------------------
Generating Weight Windows with MAGIC
------------------------------------

As discussed in the :ref:`methods section <methods_variance_reduction>`, MAGIC
is an iterative method that uses flux tally information from a Monte Carlo
simulation to produce weight windows for a user-defined mesh. While generating
the weight windows, OpenMC is capable of applying the weight windows generated
from a previous batch while processing the next batch, allowing for progressive
improvement in the weight window quality across iterations.

The typical way of generating weight windows is to define a mesh and then add an
:class:`openmc.WeightWindowGenerator` object to an :attr:`openmc.Settings`
instance, as follows::

    # Define weight window spatial mesh
    ww_mesh = openmc.RegularMesh()
    ww_mesh.dimension = (10, 10, 10)
    ww_mesh.lower_left = (0.0, 0.0, 0.0)
    ww_mesh.upper_right = (100.0, 100.0, 100.0)

    # Create weight window object and adjust parameters
    wwg = openmc.WeightWindowGenerator(
        method='magic',
        mesh=ww_mesh,
        max_realizations=settings.batches
    )

    # Add generator to Settings instance
    settings.weight_window_generators = wwg

Notably, the :attr:`max_realizations` attribute is adjusted to the number of
batches, such that all iterations are used to refine the weight window
parameters.

With the :class:`~openmc.WeightWindowGenerator` instance added to the
:attr:`~openmc.Settings`, the rest of the problem can be defined as normal. When
running, note that the second iteration and beyond may be several orders of
magnitude slower than the first. As the weight windows are applied in each
iteration, particles may be agressively split, resulting in a large number of
secondary (split) particles being generated per initial source particle. This is
not necessarily a bad thing, as the split particles are much more efficient at
exploring low flux regions of phase space as compared to initial particles.
Thus, even though the reported "particles/second" metric of OpenMC may be much
lower when generating (or just applying) weight windows as compared to analog
MC, it typically leads to an overall improvement in the figure of merit
accounting for the reduction in the variance.

.. warning::
    The number of particles per batch may need to be adjusted downward
    significantly to result in reasonable runtimes when weight windows are being
    generated or used.

At the end of the simulation, a ``weight_windows.h5`` file will be saved to disk
for later use. Loading it in another subsequent simulation will be discussed in
the "Using Weight Windows" section below.

------------------------------------------------------
Generating Weight Windows with FW-CADIS and Random Ray
------------------------------------------------------

Weight window generation with FW-CADIS and random ray in OpenMC uses the same
exact strategy as with MAGIC. An :class:`openmc.WeightWindowGenerator` object is
added to the :attr:`openmc.Settings` object, and a ``weight_windows.h5`` will be
generated at the end of the simulation. The only difference is that the code
must be run in random ray mode. A full description of how to enable and setup
random ray mode can be found in the :ref:`Random Ray User Guide <random_ray>`.

.. note::
    It is a long term goal for OpenMC to be able to generate FW-CADIS weight
    windows with only a few tweaks to an existing continuous energy Monte Carlo
    input deck. However, at the present time, the workflow requires several
    steps to generate multigroup cross section data and to configure the random
    ray solver. A high level overview of the current workflow for generation of
    weight windows with FW-CADIS using random ray is given below.

1. Produce approximate multigroup cross section data (stored in a ``mgxs.h5``
   library). There is more information on generating multigroup cross sections
   via OpenMC in the :ref:`multigroup materials <create_mgxs>` user guide, and a
   specific example of generating cross section data for use with random ray in
   the :ref:`random ray MGXS guide <mgxs_gen>`.

2. Make a copy of your continuous energy Python input file. You'll edit the new
   file to work in multigroup mode with random ray for producing weight windows.

3. Adjust the material definitions in your new multigroup Python file to utilize
   the multigroup cross sections instead of nuclide-wise continuous energy data.
   There is a specific example of making this conversion in the :ref:`random ray
   MGXS guide <mgxs_gen>`.

4. Configure OpenMC to run in random ray mode (by adding several standard random
   ray input flags and settings to the :attr:`openmc.Settings.random_ray`
   dictionary). More information can be found in the  :ref:`Random Ray User
   Guide <random_ray>`.

5. Add in a :class:`~openmc.WeightWindowGenerator` in a similar manner as for
   MAGIC generation with Monte Carlo and set the :attr:`method` attribute set to
   ``"fw_cadis"``::

       # Define weight window spatial mesh
       ww_mesh = openmc.RegularMesh()
       ww_mesh.dimension = (10, 10, 10)
       ww_mesh.lower_left = (0.0, 0.0, 0.0)
       ww_mesh.upper_right = (100.0, 100.0, 100.0)

       # Create weight window object and adjust parameters
       wwg = openmc.WeightWindowGenerator(
           method='fw_cadis',
           mesh=ww_mesh,
           max_realizations=settings.batches
       )

       # Add generator to openmc.settings object
       settings.weight_window_generators = wwg


.. warning::
    If using FW-CADIS weight window generation, ensure that the selected weight
    window mesh does not subdivide any cells in the problem. In the future, this
    restriction is intended to be relaxed, but for now subdivision of cells by a
    mesh tally will result in undefined behavior.

6. When running your multigroup random ray input deck, OpenMC will automatically
   run a forward solve followed by an adjoint solve, with a
   ``weight_windows.h5`` file generated at the end. The ``weight_windows.h5``
   file will contain FW-CADIS generated weight windows. This file can be used in
   identical manner as one generated with MAGIC, as described below.

--------------------
Using Weight Windows
--------------------

To use a ``weight_windows.h5`` weight window file with OpenMC's Monte Carlo
solver, the Python input just needs to load the h5 file::

    settings.weight_window_checkpoints = {'collision': True, 'surface': True}
    settings.survival_biasing = False
    settings.weight_windows = openmc.hdf5_to_wws('weight_windows.h5')
    settings.weight_windows_on = True

The :class:`~openmc.WeightWindowGenerator` instance is not needed to load an
existing ``weight_windows.h5`` file. Inclusion of a
:class:`~openmc.WeightWindowGenerator` instance will cause OpenMC to generate
*new* weight windows and thus overwrite any existing ``weight_windows.h5`` file.
Note that the mesh information is embedded into the weight window file, so the
mesh does not need to be redefined. Monte Carlo solves that load a weight window
file as above will utilize the weight windows to reduce the variance of the
simulation.
