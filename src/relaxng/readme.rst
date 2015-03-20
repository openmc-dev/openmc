=====================
Editing RelaxNG files
=====================

All direct edits to RelaxNG files should be in the .rnc files. The program
TRANG_ should be used to generate a correcsponding .rng file. For Ubuntu, you
can install with:

.. code-block:: bash
 
   sudo apt-get install trang

To convert the .rnc file to .rng, use the following syntax:

.. code-block:: bash
 
   trang {filename}.rnc {filename}.rng

.. _TRANG: http://www.thaiopensource.com/relaxng/trang.html
