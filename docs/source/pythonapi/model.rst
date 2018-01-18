-------------------------------------
:mod:`openmc.model` -- Model Building
-------------------------------------

Convenience Functions
---------------------

Several helper functions are available here. Ther first two create rectangular
and hexagonal prisms defined by the intersection of four and six surface
half-spaces, respectively. The last function takes a sequence of surfaces and
returns the regions that separate them.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.model.get_hexagonal_prism
   openmc.model.get_rectangular_prism
   openmc.model.subdivide

TRISO Fuel Modeling
-------------------

Classes
+++++++

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.model.TRISO

Functions
+++++++++

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.model.create_triso_lattice
   openmc.model.pack_trisos

Model Container
---------------

Classes
+++++++

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.model.Model
