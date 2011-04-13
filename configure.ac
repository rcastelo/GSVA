### configure.ac                                        -*- Autoconf -*-
###
### Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

### analogous to Perl's chop command
define([m4_chop], [m4_substr($1, 0, m4_decr(m4_len($1)))])

### get the version number from the DESCRIPTION file of the package
m4_define([GSVA_VERSION],
           m4_chop(m4_esyscmd([egrep '^Version:' DESCRIPTION | cut -d ' ' -f 2])))

AC_INIT([GSVA], [GSVA_VERSION], [justin.guinney@sagebase.org])

AC_ARG_WITH([gsl-include],
            AC_HELP_STRING([--with-gsl-include=INCLUDE_PATH],
                           [the location of the GNU Scientific Library header files]),
                           [gsl_include_path=$withval])
GSL_CPPFLAGS="-I."
if test [ -n "$gsl_include_path" ] ; then
  GSL_CPPFLAGS="-I. -I${gsl_include_path}"
else
  if test [ -n "${GSL_INCLUDE}" ] ; then
    GSL_CPPFLAGS="-I. -I${GSL_INCLUDE}"
  fi
fi

AC_ARG_WITH([gsl-lib],
            AC_HELP_STRING([--with-gsl-lib=LIB_PATH],
                           [the location of the GNU Scientific Library]),
            [gsl_lib_path=$withval])
if test [ -n "$gsl_lib_path" ] ; then
  LIBS="-L$gsl_lib_path ${LIBS}"
else
  if test [ -n "${GSL_LIBS}" ] ; then
    LIBS="-L${GSL_LIBS} ${LIBS}"
  else
    if test -n "${GSL_CONFIG}"; then
      gsl_lib_path=`gsl-config --libs | sed s/-l.*//g`
      LIBS="${gsl_lib_path} ${LIBS}"
    fi
  fi
fi

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`" -Wall -pedantic"
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC

if test -n "${GSL_CONFIG}" ; then
  GSL_CPPFLAGS=`gsl-config --cflags`
fi
CPPFLAGS="${CPPFLAGS} ${GSL_CPPFLAGS}"

if test [ -n "$gsl_lib_path" ] ; then
  LD_LIBRARY_PATH="$gsl_lib_path:${LD_LIBRARY_PATH}"
  export LD_LIBRARY_PATH
fi

AC_CHECK_LIB(gsl, gsl_stats_sd,
             [LIBS="${LIBS} -lgsl -lgslcblas"], 
             [AC_MSG_ERROR([Cannot find Gnu Scientific Library.])],
             [-lgslcblas])

AC_CHECK_HEADER(gsl/gsl_version.h, ,
             [AC_MSG_ERROR([Cannot find GSL headers.])])


AC_MSG_CHECKING([if GSL version >= 1.12])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_version.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     printf("minor=%d\n", minor);
     exit (1);
   }
   exit(minor < 12);
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 1.12])
else
  AC_MSG_RESULT([yes])
fi

AC_SUBST(GSL_CPPFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT()
