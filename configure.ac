#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(deconvolute, 1.0, helen@strubi.ox.ac.uk)
AM_INIT_AUTOMAKE(deconvolute, 1.0)
AM_INIT_AUTOMAKE(gap, 0.1)
AM_INIT_AUTOMAKE(fsmelt, 0.1)
AM_INIT_AUTOMAKE(map_to_sf, 0.1)
AC_CONFIG_SRCDIR([src/definitions.h])
AC_CONFIG_SUBDIRS([libccp4])
AC_CONFIG_HEADERS([config.h])

PACKAGE=deconvolute
VERSION=1.0

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MKDIR_P
AC_PROG_INSTALL

# Set library paths
LDFLAGS="$LDFLAGS "

# # PHENIX build/lib directory

if test -z "${LIBTBX_BUILD}" ; then
	AC_MSG_ERROR([Missing CCTBX libraries (LIBTBX_BUILD), please install CCTBX or go to the CCTBX build directory and type either: source setpaths_all.csh    (tcsh/csh) or . setpaths_all.sh (bash).])
else
	AC_MSG_RESULT([LIBTBX build directory has been found])
	CPPFLAGS="$CPPFLAGS  -I$LIBTBX_BUILD/include"
	LDFLAGS="$LDFLAGS -L$LIBTBX_BUILD/lib"
fi

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# # PHENIX ccp4io include directory

if test -z "${CCP4IO_DIST}" ; then
	AC_MSG_ERROR([Missing PHENIX libraries (CCP4IO_DIST), please install CCTBX or go to the CCTBX build directory and type either: source setpaths_all.csh    (tcsh/csh) or . setpaths_all.sh (bash).])
else
	AC_MSG_RESULT([CCP4IO directory has been found])
	CPPFLAGS="$CPPFLAGS -I$CCP4IO_DIST/libccp4/ccp4"
fi

# # PHENIX cctbx_project include directory

if test -z "${CCTBX_DIST}" ; then
	AC_MSG_ERROR([Missing PHENIX libraries (CCTBX_DIST), please install CCTBX or go to the CCTBX build directory and type either: source setpaths_all.csh    (tcsh/csh) or . setpaths_all.sh (bash).])
else
	AC_MSG_RESULT([CCTBX directory has been found])
	CPPFLAGS="$CPPFLAGS -I$CCTBX_DIST/../"
fi


# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_PROG_RANLIB
AC_PROG_F77
AC_PROG_FC

#AC_CHECK_PROG(GAP, gap, yes)
#if test x"${GAP}" != x"yes" ; then
#    AC_MSG_ERROR([Please install GAP before trying to build deconvolute.])
#fi

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt strerror])

# Check for linux 32 / 64 bit
AC_CHECK_SIZEOF([int *])

if test "${SIZEOF_INT_P}" == "4" ; then
	AC_MSG_RESULT([Configuring GAP for 32-bit operating system])
	AC_SUBST(target, linux32)
	cp gapsrc/p2-linux32.inc gapsrc/p2.inc
else
	AC_MSG_RESULT([Configuring GAP for 64-bit operating system])
	AC_SUBST(target, linux64)
	cp gapsrc/p2-linux64.inc gapsrc/p2.inc
fi

AC_OUTPUT(Makefile src/Makefile gapsrc/Makefile fsmelt/Makefile map_to_sfsrc/Makefile)
