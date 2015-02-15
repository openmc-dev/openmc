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


# AC_FC_REAL_KIND([KIND_DECLARATION], [VARIABLE_SUFFIX], [ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro checks what integer is produced by the kind 
# declaration KIND_DECLARATION. This integer is placed in 
# AC_FC_KIND_<VARIABLE_SUFFIX>. If we successfully find a
# kind integer, ACTION_IF_SUCCESS is performed; otherwise
# ACTION_IF_FAIL.

AC_DEFUN([AC_FC_REAL_KIND], [dnl
AC_REQUIRE([AC_PROG_FC])
AC_CACHE_CHECK([for kind number produced by $1], 
                ac_cv_fc_real_kind_[]$2[],
[dnl
AC_LANG_PUSH([Fortran])
FCFLAGS_save="$FCFLAGS"
FCFLAGS="$FCFLAGS $FCFLAGS_free_f90"
ac_fc_kind_test=1
ac_fc_kind_found=no
while test $ac_fc_kind_test -lt 100
do
  cat > conftest.$ac_ext << ACEOF
dnl The program below will fail to compile if 
dnl sp != mysp; ie if the kind produced by the 
dnl supplied kind declaration ($1) is not the same
dnl same as $ac_fc_kind_test. This is because Fortran
dnl pointers & targets must be of the same kind. 
dnl All conforming compilers must fail to compile the
dnl subroutine otherwise.
dnl
dnl This approach is taken since it enables us to
dnl test for kind numbers at compile time rather
dnl than run time, which means the macro will support
dnl crosss-compilation.
dnl
dnl However, note that kind numbers can theoretically
dnl be anything from 1 to the largest default integer
dnl supported by the compiler. Here we only test up to
dnl 99, which is more than enough on all compilers tried
dnl so far
dnl
      subroutine kind_explorer
      integer, parameter :: sp = $1
      integer, parameter :: mysp = $ac_fc_kind_test
      real(kind=sp), target :: x
      real(kind=mysp), pointer :: y
      y=>x
      end subroutine kind_explorer
ACEOF
dnl
  if eval $ac_compile 2>&5
  then
    ac_fc_kind_found=yes
    break
  fi
  ac_fc_kind_test=`expr $ac_fc_kind_test + 1`
done
if test "$ac_fc_kind_found" = yes; then
  ac_cv_fc_real_kind_[]$2[]=$ac_fc_kind_test
else 
  ac_cv_fc_real_kind_[]$2[]=none
fi

FCFLAGS="$FCFLAGS_save"
AC_LANG_POP([Fortran])
])
AS_IF([test $ac_cv_fc_real_kind_[]$2[] != no],
      [ac_fc_real_kind_[]$2[]=$ac_cv_fc_real_kind_[]$2[]; m4_default([$3],[:])],
      [m4_default([$4],[AC_MSG_ERROR([Could not find Fortran real kind number for $1])])]
     )
]) # AC_FC_REAL_KIND

# AC_FC_INT_KIND([KIND_DECLARATION], [VARIABLE_SUFFIX], [ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro checks what integer is produced by the kind 
# declaration KIND_DECLARATION. This integer is placed in 
# AC_FC_KIND_<VARIABLE_SUFFIX>. If we successfully find a
# kind integer, ACTION_IF_SUCCESS is performed; otherwise
# ACTION_IF_FAIL.

AC_DEFUN([AC_FC_INT_KIND], [dnl
AC_REQUIRE([AC_PROG_FC])
AC_CACHE_CHECK([for kind number produced by $1], 
                ac_cv_fc_int_kind_[]$2[],
[dnl
AC_LANG_PUSH([Fortran])
FCFLAGS_save="$FCFLAGS"
FCFLAGS="$FCFLAGS $FCFLAGS_free_f90"
ac_fc_kind_test=1
ac_fc_kind_found=no
while test $ac_fc_kind_test -lt 100
do
  cat > conftest.f90 << ACEOF
dnl The program below will fail to compile if 
dnl sp != mysp; ie if the kind produced by the 
dnl supplied kind declaration ($1) is not the same
dnl same as $ac_fc_kind_test. This is because Fortran
dnl pointers & targets must be of the same kind. 
dnl All conforming compilers must fail to compile the
dnl subroutine otherwise.
dnl
dnl This approach is taken since it enables us to
dnl test for kind numbers at compile time rather
dnl than run time, which means the macro will support
dnl crosss-compilation.
dnl
dnl However, note that kind numbers can theoretically
dnl be anything from 1 to the largest default integer
dnl supported by the compiler. Here we only test up to
dnl 99, which is more than enough on all compilers tried
dnl so far
dnl
      subroutine kind_explorer
      integer, parameter :: sp = $1
      integer, parameter :: mysp = $ac_fc_kind_test
      integer(kind=sp), target :: x
      integer(kind=mysp), pointer :: y
      y=>x
      end subroutine kind_explorer
ACEOF
dnl
  if eval $ac_compile 2>&5
  then
    ac_fc_kind_found=yes
    break
  fi
  ac_fc_kind_test=`expr $ac_fc_kind_test + 1`
done
if test "$ac_fc_kind_found" = yes; then
  ac_cv_fc_int_kind_[]$2[]=$ac_fc_kind_test
else 
  ac_cv_fc_int_kind_[]$2[]=none
fi

FCFLAGS="$FCFLAGS_save"
AC_LANG_POP([Fortran])
])
AS_IF([test $ac_cv_fc_int_kind_[]$2[] != no],
      [ac_fc_int_kind_[]$2[]=$ac_cv_fc_int_kind_[]$2[]; m4_default([$3],[:])],
      [m4_default([$4],[AC_MSG_ERROR([Could not find Fortran integer kind number for $1])])]
     )
]) # AC_FC_INT_KIND


# AC_FC_GET_REAL_KINDS([ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro attempts to find the Fortran compiler's kinds
# for the following four types of real number:
#  Compiler default (single) precision
#  Compiler double precision
#  IEEE single precision
#  IEEE double precision
# The first two are guaranteed to exist; the second two may
# or may not.
# If all 4 are succesfully found,. ACTION_IF_SUCCESS is
# performed.
# Otherwise, ACTION_IF_FAIL is performed
#
AC_DEFUN([AC_FC_GET_REAL_KINDS], [dnl
AC_REQUIRE([AC_PROG_FC])

ac_fc_got_kinds=yes

AC_FC_REAL_KIND([[kind(1.0)]], [sp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[kind(1.0d0)]], [dp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[selected_real_kind(6,34)]], [ieee_sp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[selected_real_kind(15,300)]], [ieee_dp], 
                [], [ac_got_kinds=no])

AS_IF([test $ac_fc_got_kinds != no],
      [m4_default([$1],[:])],
      [m4_default([$2],[AC_MSG_ERROR([Could not find all Fortran real kinds])])]
      )
]) dnl AC_FC_GET_REAL_KINDS

