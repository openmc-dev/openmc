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

General Rules
-------------

Conform to the Fortran 2008 standard.

Make sure code can be compiled with most common compilers, especially gfortran
and the Intel Fortran compiler. This supercedes the previous rule --- if a
Fortran 2003/2008 feature is not implemented in a common compiler, do not use
it.

Do not use special extensions that can be only be used from certain compilers.

In general, write your code in lower-case. Having code in all caps does not
enhance code readability or otherwise.

Always include comments to describe what your code is doing. Do not be afraid of
using copious amounts of comments.

Use <, >, <=, >=, ==, and /= rather than .lt., .gt., .le., .ge., .eq., and .ne.

Try to keep code within 80 columns when possible.

Don't use ``print *`` or ``write(*,*)``. If writing to a file, use a specific
unit. Writing to standard output or standard error should be handled by the
``write_message`` subroutine or functionality in the error module.

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

Variable names should be all lower-case and descriptive, i.e. not a random
assortment of letters that doesn't give any information to someone seeing it for
the first time. Variables consisting of multiple words should be separated by
underscores, not hyphens or in camel case.

Constant (parameter) variables should be in ALL CAPITAL LETTERS and defined in
in the constants.F90 module.

32-bit reals (real(4)) should never be used. Always use 64-bit reals (real(8)).

For arbitrary length character variables, use the pre-defined lengths
MAX_LINE_LEN, MAX_WORD_LEN, and MAX_FILE_LEN if possible.

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
+++++++++++++++++++++++

Always put shared variables in modules. Access module variables through a
``use`` statement. Always use the ``only`` specifier on the ``use`` statement
except for variables from the global, constants, and various header modules.

Never use ``equivalence`` statements, ``common`` blocks, or ``data`` statements.

Derived Types and Classes
-------------------------

Derived types and classes should have CamelCase names with words not separated
by underscores or hyphens.

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

------
Python
------

Style for Python code should follow PEP8_.

Docstrings for functions and methods should follow numpydoc_ style.

Python code should work with both Python 2.7+ and Python 3.0+.

Use of third-party Python packages should be limited to numpy_, scipy_, and
h5py_. Use of other third-party packages must be implemented as optional
dependencies rather than required dependencies.

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _numpydoc: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
.. _numpy: http://www.numpy.org/
.. _scipy: http://www.scipy.org/
.. _h5py: http://www.h5py.org/
