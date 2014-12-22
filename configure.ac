#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([lieeathome], [1.25])
AC_CONFIG_SRCDIR([liee/liee_worker.cpp])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

# boinc-src
AC_ARG_WITH(boinc, [  --with-boinc            enable boinc])
if test "$with_boinc" == "yes"
then
    AC_MSG_ERROR([You need to specify a path to the boinc sources as in: --with-boinc="/opt/boinc-src" ])
elif test "$with_boinc" != ""
then
    AC_MSG_NOTICE([Your path to boinc sources is $with_boinc])
    AC_SUBST(with_boinc)
else
    AC_MSG_ERROR([You need to specify a path to the boinc sources as in: --with-boinc="/opt/boinc-src" ])
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_CPP

# Checks for libraries.
# F.I.X.M.E. (for all AC_CHECK_LIB): Replace `main' with a function in the respective lib:
AC_CHECK_LIB([pthread], [main])
AC_CHECK_LIB([log4cxx], [main])
AC_CHECK_LIB([z], [main])
AC_CHECK_LIB([bz2], [main])
AC_CHECK_LIB([archive], [main])
AC_CHECK_LIB([tinyxml], [main])
AC_CHECK_LIB([boost_regex], [main])
AC_CHECK_LIB([boost_serialization], [main])
AC_CHECK_LIB([boost_system], [main])
AC_CHECK_LIB([boost_thread], [main])
AC_CHECK_LIB([boost_filesystem], [main])
#AC_CHECK_LIB([boinc], [main])
AC_CHECK_LIB([fftw3], [main])
#AC_CHECK_LIB([boinc_api], [main])
AC_CHECK_LIB([alglib], [main])
# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

#AC_CONFIG_FILES([liee/Makefile] lib/Makefile)
AC_OUTPUT(Makefile liee/Makefile lib/Makefile)
