.. _devguide_docbuild:

=============================
Building Sphinx Documentation
=============================

In order to build the documentation in the ``docs`` directory, you will need to
have the `Sphinx <http://openmc.readthedocs.io/en/latest/>`_ third-party Python
package. The easiest way to install Sphinx is via pip:

.. code-block:: sh

    sudo pip install sphinx

Additionally, you will also need a Sphinx extension for numbering figures. The
`Numfig <http://openmc.readthedocs.io/en/latest/>`_ package can be installed
directly with pip:

.. code-block:: sh

   sudo pip install sphinx-numfig

-----------------------------------
Building Documentation as a Webpage
-----------------------------------

To build the documentation as a webpage (what appears at
http://openmc.readthedocs.io), simply go to the ``docs`` directory and run:

.. code-block:: sh

    make html

-------------------------------
Building Documentation as a PDF
-------------------------------

To build PDF documentation, you will need to have a LaTeX distribution installed
on your computer. Once you have a LaTeX distribution installed, simply go to the
``docs`` directory and run:

.. code-block:: sh

     make latexpdf
