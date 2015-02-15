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

dnl Check how to get at the abort intrinsic.

AC_DEFUN([_TW_TRY_ABORT_BARE],
[
      PROGRAM TESTABORT
      CALL ABORT
      END PROGRAM TESTABORT
])
AC_DEFUN([_TW_TRY_ABORT_NAG],
[
      PROGRAM TESTABORT
      USE F90_UNIX_PROC, ONLY:ABORT
      CALL ABORT
      END PROGRAM TESTABORT
])
AC_DEFUN([_TW_TRY_ABORT_INTEL],
[
      PROGRAM TESTABORT
      CALL ABORT("")
      END PROGRAM TESTABORT
])
AC_DEFUN([_TW_TRY_ABORT_XLF],
[
      PROGRAM TESTABORT
      CALL ABORT_
      END PROGRAM TESTABORT
])

AC_DEFUN([TW_FC_CHECK_ABORT], [
AC_REQUIRE([AC_PROG_FC])dnl
dnl
AC_MSG_CHECKING([how to compile a call to ABORT])
dnl
dnl Try first with nothing
dnl
tw_abort_ok=no
dnl
dnl First check with one arg (this will fail if no args are necessary; testing
dnl in the opposite order will succeed when it shouldnt)
dnl
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_INTEL])],
    [tw_abort_ok=yes; tw_method="with argument";DEFS="$DEFS FC_HAVE_ABORT FC_ABORT_ARG"],
    [])
fi
dnl
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_XLF])],
    [tw_abort_ok=yes; tw_method="with underscore";DEFS="$DEFS FC_HAVE_ABORT FC_ABORT_UNDERSCORE"],
    [])
fi
dnl
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_BARE])],
    [tw_abort_ok=yes; tw_method=default;DEFS="$DEFS FC_HAVE_ABORT"],
    [])
fi
dnl
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_NAG])],
    [tw_abort_ok=yes; tw_method="with f90_unix_proc";DEFS="$DEFS FC_HAVE_ABORT"],
    [])
fi
dnl
dnl Cant get it to compile alone - need a compiler flag.
dnl Now try with -Vaxlib for intel:
dnl
if test $tw_abort_ok = no; then
   save_LDFLAGS=$LDFLAGS
   LDFLAGS="$LDFLAGS -Vaxlib"
   AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_BARE])],
    [tw_abort_ok=yes; tw_method="with -Vaxlib";DEFS="$DEFS FC_HAVE_ABORT"],
    [])
   if test $tw_abort_ok = no; then
      LDFLAGS=$save_LDFLAGS
   fi
fi
AC_MSG_RESULT([$tw_method])
dnl
AS_IF([test $tw_abort_ok = yes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Cannot compile call to ABORT ])])]
     )
dnl
])# TW_FC_CHECK_ABORT
