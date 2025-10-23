#                                               -*- Autoconf -*-
# borrowed from netdiffuseR configure.ac
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([ETAS], m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([Could not determine R_HOME.])
fi

CXX=$("${R_HOME}/bin/R" CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

AC_PROG_CXX
AC_LANG([C++])
AC_REQUIRE_CPP

# OpenMP detection
AC_ARG_ENABLE([openmp],
  [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
if test "x$enable_openmp" != "xno"; then
  AC_CHECK_LIB([gomp], [omp_get_num_threads],
    [OPENMP_CXXFLAGS="-fopenmp"],
    [OPENMP_CXXFLAGS=""],
    [-fopenmp])
else
  OPENMP_CXXFLAGS=""
fi
AC_SUBST([OPENMP_CXXFLAGS])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
