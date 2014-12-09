#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([lieeathome], [1.25])
AC_CONFIG_SRCDIR([liee/app_opti_wrap.cpp])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CPP

# Checks for libraries.
# F.I.X.M.E. (for all AC_CHECK_LIB): Replace `main' with a function in the respective lib:
AC_CHECK_LIB([alglib], [main])
AC_CHECK_LIB([archive], [main])
AC_CHECK_LIB([boinc], [main])
AC_CHECK_LIB([boinc_api], [main])
AC_CHECK_LIB([boost_filesystem], [main])
AC_CHECK_LIB([boost_regex], [main])
AC_CHECK_LIB([boost_serialization], [main])
AC_CHECK_LIB([boost_system], [main])
AC_CHECK_LIB([boost_thread], [main])
AC_CHECK_LIB([bz2], [main])
AC_CHECK_LIB([curl], [main])
AC_CHECK_LIB([fftw3], [main])
AC_CHECK_LIB([log4cxx], [main])
AC_CHECK_LIB([pthread], [main])
AC_CHECK_LIB([tinyxml], [main])
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

#AC_CONFIG_FILES([liee/Makefile] liee/lib/Makefile)
AC_OUTPUT(Makefile liee/Makefile liee/lib/Makefile)
