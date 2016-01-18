.. _devguide_docbuild:

=============================
Building Sphinx Documentation
=============================

In order to build the documentation in the ``docs`` directory, you will need to
have the Sphinx_ third-party Python package. The easiest way to install Sphinx
is via pip:

.. code-block:: sh

    sudo pip install sphinx

Additionally, you will also need a Sphinx extension for numbering figures. The
Numfig_ package can be installed directly with pip:

.. code-block:: sh

   sudo pip install sphinx-numfig

-----------------------------------
Building Documentation as a Webpage
-----------------------------------

To build the documentation as a webpage (what appears at
http://mit-crpg.github.io/openmc), simply go to the ``docs`` directory and run:

.. code-block:: sh

    make html

-------------------------------
Building Documentation as a PDF
-------------------------------

To build PDF documentation, you will need to have a LaTeX distribution installed
on your computer as well as Inkscape_, which is used to convert .svg files to
.pdf files. Inkscape can be installed in a Debian-derivative with:

.. code-block:: sh

    sudo apt-get install inkscape

One the pre-requisites are installed, simply go to the ``docs`` directory and
run:

.. code-block:: sh

     make latexpdf

.. _Sphinx: http://sphinx-doc.org
.. _sphinxcontrib-tikz: https://bitbucket.org/philexander/tikz
.. _Numfig: https://pypi.python.org/pypi/sphinx_numfig
.. _Inkscape: https://inkscape.org
