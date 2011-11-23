### configure.ac                                        -*- Autoconf -*-

AC_PREREQ(2.57)

AC_INIT([src/var.h.in])
AC_CONFIG_HEADERS([src/var.h])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
AC_PROG_CC

AC_CHECK_SIZEOF(double, 8)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)

AC_OUTPUT()
