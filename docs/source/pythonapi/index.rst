.. _pythonapi:

==========
Python API
==========

OpenMC includes a rich Python API that enables programmatic pre- and
post-processing. The easiest way to begin using the API is to take a look at the
:ref:`examples`. This assumes that you are already familiar with Python and
common third-party packages such as `NumPy <http://www.numpy.org/>`_. If you
have never used Python before, the prospect of learning a new code *and* a
programming language might sound daunting. However, you should keep in mind that
there are many substantial benefits to using the Python API, including:

- The ability to define dimensions using variables.
- Availability of standard-library modules for working with files.
- An entire ecosystem of third-party packages for scientific computing.
- Automated multi-group cross section generation (:mod:`openmc.mgxs`)
- A fully-featured nuclear data interface (:mod:`openmc.data`)
- Depletion capability (:mod:`openmc.deplete`)
- Convenience functions (e.g., a function returning a hexagonal region)
- Ability to plot individual universes as geometry is being created
- A :math:`k_\text{eff}` search function (:func:`openmc.search_for_keff`)
- Random sphere packing for generating TRISO particle locations
  (:func:`openmc.model.pack_spheres`)
- Ability to create materials based on natural elements or uranium enrichment

For those new to Python, there are many good tutorials available online. We
recommend going through the modules from `Codecademy
<https://www.codecademy.com/learn/learn-python-3>`_ and/or the `Scipy lectures
<https://scipy-lectures.github.io/>`_.

The full API documentation serves to provide more information on a given module
or class.

.. tip:: Users are strongly encouraged to use the Python API to generate input
         files and analyze results.

.. rubric:: Modules

.. toctree::
   :maxdepth: 1

   base
   model
   examples
   deplete
   mgxs
   stats
   data
   capi
   openmoc
