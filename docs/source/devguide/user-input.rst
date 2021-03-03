.. _devguide_user_input:

=========================
Making User Input Changes
=========================

Users are encouraged to use OpenMC's :ref:`pythonapi` to build XML files that
the OpenMC solver then reads during the initialization phase. Thus, to modify,
add, or remove user input options, changes must be made both within the Python
API and the C++ source that reads XML files produced by the Python API. The
following steps should be followed to make changes to user input:

1. Determine the Python class you need to change. For example, if you are adding
   a new setting, you probably want to change the :class:`openmc.Settings`
   class. If you are adding a new surface type, you would need to create a
   subclass of :class:`openmc.Surface`.

2. To add a new option, the class will need a `property attribute`_. For
   example, if you wanted to add a "fast_mode" setting, you would need two
   methods that look like:

   .. code-block:: python

      @property
      def fast_mode(self):
          ...

      @fast_mode.setter
      def fast_mode(self, fast_mode):
          ...

3. Make sure that when an instance of the class is exported to XML (usually
   through a ``export_to_xml()`` or ``to_xml_element()`` method), a new element
   is written to the appropriate file. OpenMC uses the
   :mod:`xml.etree.ElementTree` API, so refer to the documentation of that
   module for guidance on creating elements/attributes.

4. Make sure that your input can be categorized as one of the datatypes from
   `XML Schema Part 2`_ and that parsing of the data appropriately reflects
   this. For example, for a boolean_ value, true can be represented either by
   "true" or by "1".

5. Now that you're done with the Python side, you need to make modifications to
   the C++ codebase. Make appropriate changes in source files (e.g.,
   settings.cpp). You should use convenience functions defined by
   xml_interface.cpp.

6. If you've made changes in the geometry or materials, make sure they are
   written out to the statepoint or summary files and that the
   :class:`openmc.StatePoint` and :class:`openmc.Summary` classes read them in.

7. Finally, a set of `RELAX NG`_ schemas exists that enables validation of input
   files. You should modify the RELAX NG schema for the file you changed. The
   easiest way to do this is to change the `compact syntax`_ file
   (e.g. ``src/relaxng/geometry.rnc``) and then convert it to regular XML syntax
   using trang_::

       trang geometry.rnc geometry.rng

For most user input additions and changes, it is simple enough to follow a
"monkey see, monkey do" approach. When in doubt, contact your nearest OpenMC
developer or send a message to the `developers mailing list`_.


.. _property attribute: https://docs.python.org/3.6/library/functions.html#property
.. _XML Schema Part 2: http://www.w3.org/TR/xmlschema-2/
.. _boolean: http://www.w3.org/TR/xmlschema-2/#boolean
.. _RELAX NG: http://relaxng.org/
.. _compact syntax: http://relaxng.org/compact-tutorial-20030326.html
.. _trang: http://www.thaiopensource.com/relaxng/trang.html
.. _developers mailing list: https://groups.google.com/forum/?fromgroups=#!forum/openmc-dev
