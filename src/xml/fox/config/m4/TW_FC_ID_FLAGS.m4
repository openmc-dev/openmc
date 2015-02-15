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

dnl A macro to set various compiler-dependent things that can't be sensibly
dnl deduced.

AC_DEFUN([TW_FC_ID_FLAGS], [
AC_REQUIRE([TW_FC_ID])

case $FC_ID in

  Absoft)
     FFLAGS_DEBUG="-et -g -Rb -Rc -Rp -Rs"
     ;;

  Digital)
     FFLAGS_FAST=-O2
     FFLAGS_DEBUG="-g -Rabc -ei"
     ;;

  G77)
     ;;

  G95)
    FFLAGS_FAST=-O2
    FFLAGS_DEBUG="-ggdb3 -ftrace=full -fbounds-check -flogical=false -freal=nan -fpointer=invalid -finteger=-1"
    ;;

  Gfortran)
     ;;
 
  Intel)
     FFLAGS_DEBUG="-C -g -inline_debug_info"
     ;;

  Lahey)
     FFLAGS_DEBUG="--chk aesux --chkglobal -g --trace"
     FFLAGS_FAST="-O --warn --quiet --tpp --ntrace"
     ;;

  Nag)
     FFLAGS_DEBUG="-C=all -g -gline -nan"
     DEFS="$DEFS __NAG__"
     ;;
  
  Pathscale)
     ;;

  Portland)
     FFLAGS_DEBUG="-g -Mbounds"
     FFLAGS_FAST="-fast"
     DEFS="$DEFS PGF90"
     ;;

  SGI)
     FFLAGS_DEBUG="-g -O0"
     FFLAGS_FAST="-O3 -OPT:Olimit=0"
     ;;

  Sun)
     FFLAGS_DEBUG="-C -g"
     FFLAGS_FAST="-fast"
     ;;

  Xlf)
     FFLAGS_DEBUG="-g -C -qinitauto -qsave -qmaxmem=16000 -qnolm"
     FFLAGS_FAST="-O3 -qarch=auto -qtune=auto -qcache=auto -qnolm"
     ;;

esac

AC_SUBST(FFLAGS_MPI)

])
