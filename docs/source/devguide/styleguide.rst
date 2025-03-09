.. _devguide_styleguide:

======================
Style Guide for OpenMC
======================

In order to keep the OpenMC code base consistent in style, this guide specifies
a number of rules which should be adhered to when modified existing code or
adding new code in OpenMC.

---
C++
---

.. _styleguide_formatting:

Automatic Formatting
--------------------

To ensure consistent styling with little effort, this project uses `clang-format
<https://clang.llvm.org/docs/ClangFormat.html>`_. The repository contains a
``.clang-format`` file that can be used to automatically apply a consistent
format. The easiest way to use clang-format is to run
``tools/dev/install-commit-hooks.sh`` to install a post-commit hook that gets
executed each time a commit is made. Note that this script requires that you
already have clang-format installed. In addition, you may want to configure your
editor/IDE to automatically runs clang-format using the ``.clang-format`` file
whenever a file is saved. For example, `Visual Studio Code
<https://code.visualstudio.com/docs/cpp/cpp-ide#_code-formatting>`_ includes
support for running clang-format.

.. note::
    OpenMC's CI uses `clang-format` version 15. A different version of `clang-format`
    may produce different line changes and as a result fail the CI test.

Miscellaneous
-------------

Follow the `C++ Core Guidelines`_ except when they conflict with another
guideline listed here. For convenience, many important guidelines from that
list are repeated here.

Conform to the C++17 standard.

Always use C++-style comments (``//``) as opposed to C-style (``/**/``). (It
is more difficult to comment out a large section of code that uses C-style
comments.)

Do not use C-style casting. Always use the C++-style casts ``static_cast``,
``const_cast``, or ``reinterpret_cast``. (See `ES.49 <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#es49-if-you-must-use-a-cast-use-a-named-cast>`_)

Source Files
------------

Use a ``.cpp`` suffix for code files and ``.h`` for header files.

Header files should always use include guards with the following style (See
`SF.8 <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Rs-guards>`_):

.. code-block:: C++

    #ifndef OPENMC_MODULE_NAME_H
    #define OPENMC_MODULE_NAME_H

    namespace openmc {
    ...
    content
    ...
    }

    #endif // OPENMC_MODULE_NAME_H

Avoid hidden dependencies by always including a related header file first,
followed by C/C++ library includes, other library includes, and then local
includes. For example:

.. code-block:: C++

   // foo.cpp
   #include "foo.h"

   #include <cstddef>
   #include <iostream>
   #include <vector>

   #include "hdf5.h"
   #include "pugixml.hpp"

   #include "error.h"
   #include "random_lcg.h"

Naming
------

Struct and class names should be CamelCase, e.g. ``HexLattice``.

Functions (including member functions) should be lower-case with underscores,
e.g. ``get_indices``.

Local variables, global variables, and struct/class member variables should be
lower-case with underscores (e.g., ``n_cells``) except for physics symbols that
are written differently by convention (e.g., ``E`` for energy). Data members of
classes (but not structs) additionally have trailing underscores (e.g.,
``a_class_member_``).

The following conventions are used for variables with short names:

- ``d`` stands for "distance"
- ``E`` stands for "energy"
- ``p`` stands for "particle"
- ``r`` stands for "position"
- ``rx`` stands for "reaction"
- ``u`` stands for "direction"
- ``xs`` stands for "cross section"

All classes and non-member functions should be declared within the ``openmc``
namespace. Global variables must be declared in a namespace nested within the
``openmc`` namespace. The following sub-namespaces are in use:

- ``openmc::data``: Fundamental nuclear data (cross sections, multigroup data,
  decay constants, etc.)
- ``openmc::model``: Variables related to geometry, materials, and tallies
- ``openmc::settings``: Global settings / options
- ``openmc::simulation``: Variables used only during a simulation

Accessors and mutators (get and set functions) may be named like
variables. These often correspond to actual member variables, but this is not
required. For example, ``int count()`` and ``void set_count(int count)``.

Variables declared constexpr or const that have static storage duration (exist
for the duration of the program) should be upper-case with underscores,
e.g., ``SQRT_PI``.

Documentation
-------------

Classes, structs, and functions are to be annotated for the `Doxygen
<https://www.doxygen.nl/>`_ documentation generation tool. Use the ``\`` form of
Doxygen commands, e.g., ``\brief`` instead of ``@brief``.

------
Python
------

Style for Python code should follow PEP8_.

Docstrings for functions and methods should follow numpydoc_ style.

Python code should work with Python 3.8+.

Use of third-party Python packages should be limited to numpy_, scipy_,
matplotlib_, pandas_, and h5py_. Use of other third-party packages must be
implemented as optional dependencies rather than required dependencies.

Prefer pathlib_ when working with filesystem paths over functions in the os_
module or other standard-library modules. Functions that accept arguments that
represent a filesystem path should work with both strings and Path_ objects.

.. _C++ Core Guidelines: http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines
.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _numpydoc: https://numpydoc.readthedocs.io/en/latest/format.html
.. _numpy: https://numpy.org/
.. _scipy: https://www.scipy.org/
.. _matplotlib: https://matplotlib.org/
.. _pandas: https://pandas.pydata.org/
.. _h5py: https://www.h5py.org/
.. _pathlib: https://docs.python.org/3/library/pathlib.html
.. _os: https://docs.python.org/3/library/os.html
.. _Path: https://docs.python.org/3/library/pathlib.html#pathlib.Path
