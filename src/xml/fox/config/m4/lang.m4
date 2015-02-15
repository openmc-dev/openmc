# _AC_COMPILER_EXEEXT_DEFAULT
# ---------------------------
# Check for the extension used for the default name for executables.
#
# We do this in order to find out what is the extension we must add for
# creating executables (see _AC_COMPILER_EXEEXT's comments).
#
# Beware of `expr' that may return `0' or `'.  Since this macro is
# the first one in touch with the compiler, it should also check that
# it compiles properly.
#
# On OpenVMS 7.1 system, the DEC C 5.5 compiler when called through a
# GNV (gnv.sourceforge.net) cc wrapper, produces the output file named
# `a_out.exe'.
m4_define([_AC_COMPILER_EXEEXT_DEFAULT],
[# First try to determine the flag needed to name the executable
# It is nearly always "-o" but Lahey Fortran wants "-out"
AC_MSG_CHECKING([for linker flag to name executables])
for ac_link_obj_flag in "/exe:" "-out " "-o "; do
AS_IF([_AC_DO_VAR(ac_link)], 
[ac_link_obj_flag_found=yes; break], 
[:])
done
if test x$ac_link_obj_flag_found = x ; then
AC_MSG_FAILURE([Could not determine flag to name executables])
fi
AC_MSG_RESULT([$ac_link_obj_flag])

# Try to create an executable without -o first, disregard a.out.
# It will help us diagnose broken compilers, and finding out an intuition
# of exeext.
AC_MSG_CHECKING([for _AC_LANG compiler default output file name])
ac_link_default=`echo "$ac_link" | sed ['s/ $ac_link_obj_flag *conftest[^ ]*//']`
#
# List of possible output files, starting from the most likely.
# The algorithm is not robust to junk in `.', hence go to wildcards (a.*)
# only as a last resort.  b.out is created by i960 compilers.
ac_files='a_out.exe a.exe conftest.exe a.out conftest a.* conftest.* b.out'
#
# The IRIX 6 linker writes into existing files which may not be
# executable, retaining their permissions.  Remove them first so a
# subsequent execution test works.
ac_rmfiles=
for ac_file in $ac_files
do
  case $ac_file in
    _AC_COMPILER_EXEEXT_REJECT ) ;;
    * ) ac_rmfiles="$ac_rmfiles $ac_file";;
  esac
done
rm -f $ac_rmfiles

AS_IF([_AC_DO_VAR(ac_link_default)],
[# Autoconf-2.13 could set the ac_cv_exeext variable to `no'.
# So ignore a value of `no', otherwise this would lead to `EXEEXT = no'
# in a Makefile.  We should not override ac_cv_exeext if it was cached,
# so that the user can short-circuit this test for compilers unknown to
# Autoconf.
for ac_file in $ac_files
do
  test -f "$ac_file" || continue
  case $ac_file in
    _AC_COMPILER_EXEEXT_REJECT )
        ;; 
    [[ab]].out )
        # We found the default executable, but exeext='' is most
        # certainly right.
        break;;
    *.* )
        if test "${ac_cv_exeext+set}" = set && test "$ac_cv_exeext" != no;
        then :; else
           ac_cv_exeext=`expr "$ac_file" : ['[^.]*\(\..*\)']`
        fi
        # We set ac_cv_exeext here because the later test for it is not
        # safe: cross compilers may not add the suffix if given an `-o'
        # argument, so we may need to know it at that point already.
        # Even if this section looks crufty: it has the advantage of
        # actually working.
        break;;
    * )
        break;; 
  esac 
done
test "$ac_cv_exeext" = no && ac_cv_exeext=
],
      [_AC_MSG_LOG_CONFTEST
AC_MSG_FAILURE([_AC_LANG compiler cannot create executables], 77)])
ac_exeext=$ac_cv_exeext
AC_MSG_RESULT([$ac_file])
])# _AC_COMPILER_EXEEXT_DEFAULT
