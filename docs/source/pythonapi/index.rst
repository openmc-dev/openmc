.. _pythonapi:

==========
Python API
==========

OpenMC includes a rich Python API that enables programmatic pre- and
post-processing. The easiest way to begin using the API is to take a look at the
example Jupyter_ notebooks provided. However, this assumes that you are already
familiar with Python and common third-party packages such as NumPy_. If you have
never programmed in Python before, there are many good tutorials available
online. We recommend going through the modules from Codecademy_ and/or the
`Scipy lectures`_. The full API documentation serves to provide more information
on a given module or class.

**Handling nuclear data:**

.. toctree::
    :maxdepth: 1

    ace
    mgxs_library

**Creating input files:**

.. toctree::
    :maxdepth: 1

    cmfd
    element
    filter
    geometry
    material
    mesh
    nuclide
    opencg_compatible
    plots
    settings
    source
    stats
    surface
    tallies
    trigger
    universe

**Running OpenMC:**

.. toctree::
    :maxdepth: 1

    executor

**Post-processing:**

.. toctree::
    :maxdepth: 1

    particle_restart
    statepoint
    summary
    tallies

**Multi-Group Cross Section Generation**

.. toctree::
    :maxdepth: 1

    mgxs

**Example Jupyter Notebooks:**

.. toctree::
    :maxdepth: 1

    examples/post-processing
    examples/pandas-dataframes
    examples/tally-arithmetic
    examples/mgxs-part-i
    examples/mgxs-part-ii
    examples/mgxs-part-iii

.. _Jupyter: https://jupyter.org/
.. _NumPy: http://www.numpy.org/
.. _Codecademy: https://www.codecademy.com/tracks/python
.. _Scipy lectures: https://scipy-lectures.github.io/
