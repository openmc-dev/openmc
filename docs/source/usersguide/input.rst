.. _usersguide_input:

========================
Creating XML Input Files
========================

Unlike many other Monte Carlo codes which use an arbitrary-format ASCII file
with "cards" to specify a particular geometry, materials, and associated run
settings, the input files for OpenMC are structured in a set of XML_ files. XML,
which stands for eXtensible Markup Language, is a simple format that allows data
to be exchanged efficiently between different programs and interfaces.

Anyone who has ever seen webpages written in HTML will be familiar with the
structure of XML whereby "tags" enclosed in angle brackets denote that a
particular piece of data will follow. Let us examine the follow example::

    <person>
      <firstname>John</firstname>
      <lastname>Smith</lastname>
      <age>27</age>
      <occupation>Health Physicist</occupation>
    </person>

Here we see that the first tag indicates that the following data will describe a
person. The nested tags *firstname*, *lastname*, *age*, and *occupation*
indicate characteristics about the person being described.

In much the same way, OpenMC input uses XML tags to describe the geometry, the
materials, and settings for a Monte Carlo simulation.

.. _XML: http://www.w3.org/XML/

-----
Files
-----

To assemble a complete model for OpenMC, one needs to create separate XML files
for the geometry, materails, and settings. Additionally, an optional tallies XML
file specifies physical quantities to be tallied. OpenMC expects that these
files are called:

* ``geometry.xml``
* ``materials.xml``
* ``setings.xml``
* ``tallies.xml``

--------------------------------------
Geometry Specification -- geometry.xml
--------------------------------------

Types of surfaces:

x-plane
  A plane perpendicular to the x axis, i.e. a surface of the form x - x0
  = 0. The coefficients specified are "x0".

y-plane
  A plane perpendicular to the y axis, i.e. a surface of the form y - y0
  = 0. The coefficients specified are "y0".

z-plane
  A plane perpendicular to the z axis, i.e. a surface of the form z - z0
  = 0. The coefficients specified are "z0".

plane
  An arbitrary plane of the form A*x + B*y + C*z = D. The coefficients
  specified are "A B C D".

x-cylinder
  An infinite cylinder whose length is paralle to the x-axis. This is a
  quadratic surface of the form (y - y0)^2 + (z - z0)^2 = R^2. The coefficients
  specified are "y0 z0 R".

y-cylinder
  An infinite cylinder whose length is paralle to the y-axis. This is a
  quadratic surface of the form (x - x0)^2 + (z - z0)^2 = R^2. The coefficients
  specified are "x0 z0 R".

z-cylinder
  An infinite cylinder whose length is paralle to the z-axis. This is a
  quadratic surface of the form (x - x0)^2 + (y - y0)^2 = R^2. The coefficients
  specified are "x0 y0 R".

sphere 
  A sphere of the form (x - x0)^2 + (y - y0)^2 + (z - z0)^2 = R^2. The
  coefficients specified are "x0 y0 z0 R".

----------------------------------------
Materials Specification -- materials.xml
----------------------------------------

--------------------------------------
Settings Specification -- settings.xml
--------------------------------------

------------------------------------
Tallies Specification -- tallies.xml
------------------------------------

