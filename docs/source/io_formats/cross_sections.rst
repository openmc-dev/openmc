.. _io_cross_sections:

============================================
Cross Sections Listing -- cross_sections.xml
============================================

.. _directory_element:

-----------------------
``<directory>`` Element
-----------------------

The ``<directory>`` element specifies a root directory to which the path for all
files listed in a :ref:`library_element` are given relative to. This element has
no attributes or sub-elements; the directory should be given within the text
node. For example,

.. code-block:: xml

   <directory>/opt/data/cross_sections/</directory>

.. _library_element:

---------------------
``<library>`` Element
---------------------

The ``<library>`` element indicates where an HDF5 cross section file is located,
whether it contains incident neutron or thermal scattering data, and what
materials are listed within. It has the following attributes:

  :materials:

    A space-separated list of nuclides or thermal scattering tables. For
    example,

    .. code-block:: xml

       <library materials="U234 U235 U238" />
       <library materials="c_H_in_H2O c_D_in_G2O" />

    Often, just a single nuclide or thermal scattering table is contained in a
    given file.

  :path:
    Path to the HDF5 file. If the :ref:`directory_element` is specified, the
    path is relative to the directory given. Otherwise, it is relative to the
    directory containing the ``cross_sections.xml`` file.

  :type:
    The type of data contained in the file, either 'neutron' or 'thermal'.
