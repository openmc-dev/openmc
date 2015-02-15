dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006

dnl @synopsis TW_CHECK_FC_95([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with Fortran 95 (ISO-IEC 1539-1:1997)
dnl 
dnl It currently tests for:
dnl
dnl Named End Interface
dnl Derived type initialization
dnl The Null() intrinsic
dnl The Forall statement 
dnl The Cpu_Time intrinsic
dnl Pure functions
dnl Elemental functions
dnl 
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_95],[
dnl
AC_MSG_CHECKING([for Fortran 95 compliance])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_f95

!      Interface test_interface
!      End Interface test_interface

      Type test_type
        Integer :: i = 1
      End Type test_type

      Integer, Pointer :: j => Null()

      Integer :: i
      Real :: a

      Forall (i=1:50)
      End Forall

      Call CPU_TIME(a)

      Contains

      Pure Integer Function test_pure()
        test_pure = 0
      End Function test_pure

      Elemental Integer Function test_elemental(in)
        Integer, Intent(In) :: in
        test_elemental = 0
      End Function test_elemental

      End Program test_f95
   ]]),
   [AC_MSG_RESULT([yes])
    m4_default([$1],[:])
   ],
   [AC_MSG_RESULT([no])
    m4_default([$2], 
               [AC_MSG_ERROR([A fully Fortran-95-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
