==========================================
OpenMC Monte Carlo Particle Transport Code
==========================================

|licensebadge| |travisbadge| |coverallsbadge|

The OpenMC project aims to provide a fully-featured Monte Carlo particle
transport code based on modern methods. It is a constructive solid geometry,
continuous-energy transport code that uses HDF5 format cross sections. The
project started under the Computational Reactor Physics Group at MIT.

Complete documentation on the usage of OpenMC is hosted on Read the Docs (both
for the `latest release`_ and developmental_ version). If you are interested in
the project or would like to help and contribute, please send a message to the
OpenMC User's Group `mailing list`_.

------------
Installation
------------

Detailed `installation instructions`_ can be found in the User's Guide.

------
Citing
------

If you use OpenMC in your research, please consider giving proper attribution by
citing the following publication:

- Paul K. Romano, Nicholas E. Horelik, Bryan R. Herman, Adam G. Nelson, Benoit
  Forget, and Kord Smith, "`OpenMC: A State-of-the-Art Monte Carlo Code for
  Research and Development <https://doi.org/10.1016/j.anucene.2014.07.048>`_,"
  *Ann. Nucl. Energy*, **82**, 90--97 (2015).

---------------
Troubleshooting
---------------

If you run into problems compiling, installing, or running OpenMC, first check
the `Troubleshooting section`_ in the User's Guide. If you are not able to find
a solution to your problem there, please send a message to the User's Group
`mailing list`_.

--------------
Reporting Bugs
--------------

OpenMC is hosted on GitHub and all bugs are reported and tracked through the
Issues_ feature on GitHub. However, GitHub Issues should not be used for common
troubleshooting purposes. If you are having trouble installing the code or
getting your model to run properly, you should first send a message to the
User's Group `mailing list`_. If it turns out your issue really is a bug in the
code, an issue will then be created on GitHub. If you want to request that a
feature be added to the code, you may create an Issue on github.

-------
License
-------

OpenMC is distributed under the MIT/X license_.

.. _latest release: http://openmc.readthedocs.io/en/stable/
.. _developmental: http://openmc.readthedocs.io/en/latest/
.. _mailing list: https://groups.google.com/forum/?fromgroups=#!forum/openmc-users
.. _installation instructions: http://openmc.readthedocs.io/en/stable/usersguide/install.html
.. _Troubleshooting section: http://openmc.readthedocs.io/en/stable/usersguide/troubleshoot.html
.. _Issues: https://github.com/mit-crpg/openmc/issues
.. _license: http://openmc.readthedocs.io/en/stable/license.html

.. |licensebadge| image:: https://img.shields.io/github/license/mit-crpg/openmc.svg
   :target: http://openmc.readthedocs.io/en/latest/license.html
   :alt: License

.. |travisbadge| image:: https://travis-ci.org/mit-crpg/openmc.svg?branch=develop
   :target: https://travis-ci.org/mit-crpg/openmc
   :alt: Travis CI build status (Linux)

.. |coverallsbadge| image:: https://coveralls.io/repos/github/mit-crpg/openmc/badge.svg?branch=develop
   :target: https://coveralls.io/github/mit-crpg/openmc?branch=develop
   :alt: Code Coverage
