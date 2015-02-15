.SUFFIXES:
.SUFFIXES: ac m4

M4_FILES=m4/TW_CHECK_FC_90.m4 m4/TW_CHECK_FC_95.m4 m4/TW_CHECK_FC_FPP.m4 \
         m4/TW_CHECK_FC_FPP_90.m4 m4/TW_CHECK_FC_TR15580.m4 \
         m4/TW_CHECK_FC_TR15581.m4 m4/TW_FC_CHECK_ABORT.m4 \
         m4/TW_FC_CHECK_FLUSH.m4 m4/TW_FC_ID.m4 m4/TW_FC_ID_FLAGS.m4 \
         m4/TW_FC_KINDS.m4 m4/TW_PROG_CYGPATH_W.m4 m4/TW_FC_CHECK_EOL.m4 \
         m4/fortran.m4 m4/TW_PATH_NETCDF.m4 m4/lang.m4 \
         m4/TW_CHECK_ASSOCIATED_BUG.m4

../configure: configure.ac aclocal.m4
	autoconf
	mv configure ..

aclocal.m4: $(M4_FILES)
	cat m4/*.m4 > aclocal.m4
