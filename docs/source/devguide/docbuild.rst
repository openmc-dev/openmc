.. _devguide_docbuild:

=============================
Building Sphinx Documentation
=============================

In order to build the documentation in the ``docs`` directory, you will need to
have the several third-party Python packages installed, including `Sphinx
<https://www.sphinx-doc.org/en/master/>`_. To install the necessary
prerequisites, provide the optional "docs" dependencies when installing OpenMC's
Python API. That is, from the root directory of the OpenMC repository:

.. code-block:: sh

    python -m pip install ".[docs]"

The OpenMC documentation also uses Doxygen to automatically generate its
C/C++ API documentation directly from the docstrings available in the source
code. You will need to have a working installation of Doxygen to generate the
documentation locally.

-----------------------------------
Building Documentation as a Webpage
-----------------------------------

To build the documentation as a webpage (what appears at
https://docs.openmc.org), simply go to the ``docs`` directory and run:

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
