This directory contains the data needed in order to run the tests using the NDPP_ tallies.

The files in this directory are as follows, with a brief description provided:
- cross_sections.xml: A version of the cross_sections.XML file generated with the NNDC_ script in the openmc/data/ directory.
This file is used by NDPP to determine which nuclides (and temperatures) to process.
In the case of these tests, only room temperature H-1, O-16, U-235, and the H-H2O S(a,b) table are considered.
This file is included merely so the NDPP results can be re-generated in the future as needs dictate.

- ndpp.xml: The options given to NDPP to generate this data.  
Note that the options are chosen to have a small NDPP run time, not to provide accurate results.
This file is included merely so the NDPP results can be re-generated in the future as needs dictate.

- ndpp_lib.xml: The file describing the data created by NDPP so that it can be utilized by OpenMC.
This file is an output of NDPP.

- 1001.71c.g2: The 2-group data from NDPP for H-1 at room temperature.
This file is in ASCII format.

- 8016.71c.g2: The 2-group data from NDPP for O-16 at room temperature.
This file is in ASCII format.

- 92235.71c.g2: The 2-group data from NDPP for U-235 at room temperature.
This file is in ASCII format.

- HH20.71c.g2: The 2-group data from NDPP for the H-H2O S(a,b) table at room temperature.
This file is in ASCII format.

.. _NDPP:  http://ndpp.github.io/
.. _NNDC:  http://www.nndc.bnl.gov/endf/b7.1/acefiles.html
