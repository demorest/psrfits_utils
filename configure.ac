AC_INIT([psrfits_utils], [1.0])
AC_CONFIG_SRCDIR([downsample.c])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CC
AC_PROG_LIBTOOL

AX_PTHREAD

#AC_HAVE_LIBRARY([cfitsio])
AC_CHECK_LIB([cfitsio],[ffgky])

#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
	Makefile
])

AC_OUTPUT
