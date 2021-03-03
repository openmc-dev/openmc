.. _devguide_docbuild:

=============================
Building Sphinx Documentation
=============================

In order to build the documentation in the ``docs`` directory, you will need to
have the `Sphinx <https://www.sphinx-doc.org/en/master/>`_ third-party Python
package. The easiest way to install Sphinx is via pip:

.. code-block:: sh

    pip install sphinx

Additionally, you will need several Sphinx extensions that can be installed
directly with pip:

.. code-block:: sh

    pip install sphinx-numfig
    pip install sphinxcontrib-katex
    pip install sphinxcontrib-svg2pdfconverter

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
