.. _devguide_styleguide:

======================
Style Guide for OpenMC
======================

In order to keep the OpenMC code base consistent in style, this guide specifies
a number of rules which should be adhered to when modified existing code or
adding new code in OpenMC.

-------
Fortran
-------

Miscellaneous
-------------

Conform to the Fortran 2008 standard.

Make sure code can be compiled with most common compilers, especially gfortran
and the Intel Fortran compiler. This supersedes the previous rule --- if a
Fortran 2003/2008 feature is not implemented in a common compiler, do not use
it.

Do not use special extensions that can be only be used from certain compilers.

Always include comments to describe what your code is doing. Do not be afraid of
using copious amounts of comments.

Use <, >, <=, >=, ==, and /= rather than .lt., .gt., .le., .ge., .eq., and .ne.

Try to keep code within 80 columns when possible.

Don't use ``print *`` or ``write(*,*)``. If writing to a file, use a specific
unit. Writing to standard output or standard error should be handled by the
``write_message`` subroutine or functionality in the error module.

Naming
------

In general, write your code in lower-case. Having code in all caps does not
enhance code readability or otherwise.

Module names should be lower-case with underscores if needed, e.g.
``xml_interface``.

Class names should be CamelCase, e.g. ``HexLattice``.

Functions and subroutines (including type-bound methods) should be lower-case
with underscores, e.g. ``get_indices``.

Local variables, global variables, and type attributes should be lower-case
with underscores (e.g. ``n_cells``) except for physics symbols that are written
differently by convention (e.g. ``E`` for energy).

Constant (parameter) variables should be in upper-case with underscores, e.g.
``SQRT_PI``. If they are used by more than one module, define them in the
constants.F90 module.

Procedures
----------

Above each procedure, include a comment block giving a brief description of what
the procedure does.

Nonpointer dummy arguments to procedures should be explicitly specified as
intent(in), intent(out), or intent(inout).

Include a comment describing what each argument to a procedure is.

Variables
---------

Never, under any circumstances, should implicit variables be used! Always
include ``implicit none`` and define all your variables.

32-bit reals (real(4)) should never be used. Always use 64-bit reals (real(8)).

For arbitrary length character variables, use the pre-defined lengths
``MAX_LINE_LEN``, ``MAX_WORD_LEN``, and ``MAX_FILE_LEN`` if possible.

Do not use old-style character/array length (e.g. character*80, real*8).

Integer values being used to indicate a certain state should be defined as named
constants (see the constants.F90 module for many examples).

Always use a double colon :: when declaring a variable.

Yes:

.. code-block:: fortran

    if (boundary_condition == BC_VACUUM) then

No:

.. code-block:: fortran

    if (boundary_condition == -10) then

Avoid creating arrays with a pre-defined maximum length. Use dynamic memory
allocation instead. Use allocatable variables instead of pointer variables when
possible.

Shared/Module Variables
-----------------------

Always put shared variables in modules. Access module variables through a
``use`` statement. Always use the ``only`` specifier on the ``use`` statement
except for variables from the global, constants, and various header modules.

Never use ``equivalence`` statements, ``common`` blocks, or ``data`` statements.

Indentation
-----------

Never use tab characters. Indentation should always be applied using
spaces. Emacs users should include the following line in their .emacs file:

.. code-block:: common-lisp

    (setq-default indent-tabs-mode nil)

vim users should include the following line in their .vimrc file::

    set expandtab

Use 2 spaces per indentation level. This applies to all constructs such as
program, subroutine, function, if, associate, etc. Emacs users should set the
variables f90-if-indent, f90-do-indent, f90-continuation-indent,
f90-type-indent, f90-associate-indent, and f90-program indent to 2.

Continuation lines should be indented by at least 5 spaces. They may be indented
more in order to make the content match the context.  For example, either of
these are valid continuation indentations:

.. code-block:: fortran

    local_xyz(1) = xyz(1) - (this % lower_left(1) + &
         (i_xyz(1) - HALF)*this % pitch(1))
    call which_data(scatt_type, get_scatt, get_nuscatt, get_chi_t, get_chi_p, &
                    get_chi_d, scatt_order)

Whitespace in Expressions
-------------------------

Use a single space between arguments to procedures.

Avoid extraneous whitespace in the following situations:

- In procedure calls::

    Yes: call somesub(x, y(2), z)
    No:  call somesub( x, y( 2 ), z )

- In logical expressions, use one space around operators but nowhere else::

    Yes: if (variable == 2) then
    No:  if ( variable==2 ) then

The structure component designator ``%`` should be surrounded by one space on
each side.

Do not leave trailing whitespace at the end of a line.

---
C++
---

Miscellaneous
-------------

Follow the `C++ Core Guidelines`_ except when they conflict with another
guideline listed here. For convenience, many important guidelines from that
list are repeated here.

Conform to the C++11 standard. Note that this is a significant difference
between our style and the C++ Core Guidelines.  Many suggestions in those
Guidelines require C++14.

Always use C++-style comments (``//``) as opposed to C-style (``/**/``). (It
is more difficult to comment out a large section of code that uses C-style
comments.)

Header files should always use include guards with the following style (See
`SF.8 <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#sf8-use-include-guards-for-all-h-files>`_:

.. code-block:: C++

    #ifndef MODULE_NAME_H
    #define MODULE_NAME_H
    ...
    content
    ...
    #endif // MODULE_NAME_H

Do not use C-style casting. Always use the C++-style casts ``static_cast``,
``const_cast``, or ``reinterpret_cast``. (See `ES.49 <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#es49-if-you-must-use-a-cast-use-a-named-cast>`_)

Naming
------

In general, write your code in lower-case. Having code in all caps does not
enhance code readability or otherwise.

Struct and class names should be CamelCase, e.g. ``HexLattice``.

Functions (including member functions) should be lower-case with underscores,
e.g. ``get_indices``.

Local variables, global variables, and struct/class attributes should be
lower-case with underscores (e.g. ``n_cells``) except for physics symbols that
are written differently by convention (e.g. ``E`` for energy).

Const variables should be in upper-case with underscores, e.g. ``SQRT_PI``.

Curly braces
------------

For a function definition, the opening and closing braces should each be on
their own lines.  This helps distinguish function code from the argument list.
If the entire function fits on one line, then the braces can be on the same
line. e.g.:

.. code-block:: C++

    return_type function(type1 arg1, type2 arg2)
    {
      content();
    }

    return_type
    function_with_many_args(type1 arg1, type2 arg2, type3 arg3,
                            type4 arg4)
    {
      content();
    }

    int return_one() {return 1;}

For a conditional, the opening brace should be on the same line as the end of
the conditional statement. If there is a following ``else if`` or ``else``
statement, the closing brace should be on the same line as that following
statement. Otherwise, the closing brace should be on its own line. A one-line
conditional can have the closing brace on the same line or it can omit the
braces entirely e.g.:

.. code-block:: C++

    if (condition) {
      content();
    }

    if (condition1) {
      content();
    } else if (condition 2) {
      more_content();
    } else {
      further_content();
    }

    if (condition) {content()};

    if (condition) content();

For loops similarly have an opening brace on the same line as the statement and
a closing brace on its own line. One-line loops may have the closing brace on
the same line or omit the braces entirely.

.. code-block:: C++

    for (int i = 0; i < 5; i++) {
      content();
    }

    for (int i = 0; i < 5; i++) {content();}

    for (int i = 0; i < 5; i++) content();

------
Python
------

Style for Python code should follow PEP8_.

Docstrings for functions and methods should follow numpydoc_ style.

Python code should work with both Python 2.7+ and Python 3.0+.

Use of third-party Python packages should be limited to numpy_, scipy_, and
h5py_. Use of other third-party packages must be implemented as optional
dependencies rather than required dependencies.

.. _C++ Core Guidelines: http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines
.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _numpydoc: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
.. _numpy: http://www.numpy.org/
.. _scipy: http://www.scipy.org/
.. _h5py: http://www.h5py.org/
