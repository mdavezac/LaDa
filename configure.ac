AC_INIT(lamarck/structure.cc)
AM_INIT_AUTOMAKE(lamarck,0.1)
AC_PROG_CC
AC_PROG_CXX
#AC_PROG_F77(g95)
#AC_PROG_FC(g95)
AC_PROG_INSTALL

AM_CONFIG_HEADER(config.h)
# TinyXml shit
AH_TOP([
#ifndef _G_HAVE_BOOL
#define _G_HAVE_BOOL
#endif])
AH_TOP([
#ifndef ANSI_HEADERS
#define ANSI_HEADERS
#endif])
AH_TOP([
#ifndef HAVE_SSTREAM
#define HAVE_SSTREAM
#endif])


AC_ARG_ENABLE(vff,
     [  --enable-vff    Turn on vff],
     [case "${enableval}" in
       yes) vff=true ;;
       no)  vff=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-vff) ;;
     esac],[vff=true])

AC_ARG_ENABLE(pescan,
     [  --enable-pescan    Turn on pescan],
     [case "${enableval}" in
       yes) pescan=true ;;
       no)  pescan=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-pescan) ;;
     esac],[pescan=true])

AC_ARG_ENABLE(darwin,
     [  --enable-darwin    Turn on darwin],
     [case "${enableval}" in
       yes) darwin=true ;;
       no)  darwin=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-darwin) ;;
     esac],[darwin=true])

AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=true])
if test x$debug = xtrue; then
       CPP_DEBUG_FLAGS="-Wall -O0  -fbounds-check"
     else
       CPP_DEBUG_FLAGS=
fi
AC_SUBST(CPP_DEBUG_FLAGS)

AC_ARG_ENABLE(mpi,
[  --enable-mpi    Turn on mpi],
[case "${enableval}" in
  yes) mpi=true ;;
  no)  mpi=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
esac],[mpi=false])
AH_TEMPLATE([_MPI])

#  MPI libs
if test x$mpi = xtrue; then
  PKG_CHECK_MODULES(MPICH, mpich >= 1.2.7 )
  AC_DEFINE([_MPI],[1],[])
  MPI_ADDLIB="../mpi/libmpiobject.la"
fi
AC_SUBST(MPICH_CFLAGS)
AC_SUBST(MPICH_LIBS)
AC_SUBST(MPI_ADDLIB)

#  XML libs
PKG_CHECK_MODULES(TINYXML, tinyxml)
AC_SUBST(TINYXML_LIBS)
AC_SUBST(TINYXML_CFLAGS)

#  libm 
AC_CHECK_LIB(m,main)

# enables vff if pescan and darwin are enabled is true
if test x$pescan = xtrue && test x$darwin = xtrue; then
  $vff = true
fi

# conditional
AM_CONDITIONAL(DOMPI, test x$mpi = xtrue)
AM_CONDITIONAL(DOPESCAN, test x$pescan = xtrue)
AM_CONDITIONAL(DOVFF, test x$vff = xtrue)
AM_CONDITIONAL(DODARWIN, test x$darwin = xtrue)

#  gsl libs iff vff is enabled 
if test x$vff = xtrue; then
  AC_FIND_FILE(libgslcblas.a,[$HOME/usr/lib /usr/lib /usr/local/lib], gsllibcblas)
  AC_FIND_FILE(libgsl.a,[$HOME/usr/lib /usr/lib /usr/local/lib], gsllib)
  AC_FIND_FILE(gsl_multimin.h,[$HOME/usr/include/gsl /usr/include/gsl /usr/local/include/gsl], gslinclude)
  if test $gslinclude != NO; then
    GSL_INCLUDE="-I $gslinclude"
    GSL_LIB="-L$gsllibcblas -lgslcblas -L$gsllib -lgsl"
  else
    AC_MSG_ERROR("Could not find gsl headers")
  fi
  if test $gsllib = NO; then
    AC_MSG_ERROR("Could not find libgsl.a")
  fi
  if test $gsllibcblas = NO; then
    AC_MSG_ERROR("Could not find libgslcblas.a")
  fi
fi 
AC_SUBST(GSL_INCLUDE)
AC_SUBST(GSL_LIB)

#  eo libs iff darwin is enabled 
if test x$darwin = xtrue; then
  PKG_CHECK_MODULES(EO, eo >= 1.0.1 )
fi 
AC_SUBST(GSL_INCLUDE)
AC_SUBST(GSL_LIB)

# checks for STLs
AC_CXX_HAVE_STL
if test "$ac_cv_cxx_have_stl" = yes; then
  LIBS="-lstdc++ -lm $LIBS"
 else
   AC_MSG_ERROR("Could not find standard template libraries")
fi

# checks for compiler specific flags 
#
COMPILER_CXXFLAGS=
# first checks if friend injection is needed (for atat) (gcc >= 4.0)
if test "$ac_compiler_gnu"="yes"; then
AX_CHECK_COMPILER_FLAGS([-ffriend-injection], 
                        [COMPILER_CXXFLAGS+=" -ffriend-injection"],)
AX_CHECK_COMPILER_FLAGS([-malign-double], 
                        [COMPILER_CXXFLAGS+=" -malign-double"],)
fi
AC_SUBST(COMPILER_CXXFLAGS)

#
AC_PROG_LIBTOOL


AC_OUTPUT(Makefile atat/Makefile physics/Makefile opt/Makefile lamarck/Makefile
          vff/Makefile mpi/Makefile pescan_interface/Makefile darwin/Makefile )


