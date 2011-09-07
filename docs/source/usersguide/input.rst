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

----------------------------------------
Materials Specification -- materials.xml
----------------------------------------

--------------------------------------
Settings Specification -- settings.xml
--------------------------------------

------------------------------------
Tallies Specification -- tallies.xml
------------------------------------

