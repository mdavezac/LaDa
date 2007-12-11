AC_INIT(lamarck/structure.cc)
AM_INIT_AUTOMAKE(LaDa,3.0)
# if test $F77 == ""; then
#   echo "HERE\n\n";
# fi

AC_PROG_CC
AC_PROG_CXX

AC_PROG_INSTALL

AM_CONFIG_HEADER(config.h)
AH_TOP([
#ifndef _G_HAVE_BOOL
#define _G_HAVE_BOOL
#endif
#ifndef ANSI_HEADERS
#define ANSI_HEADERS
#endif
#ifndef HAVE_SSTREAM
#define HAVE_SSTREAM
#endif])

AC_ARG_ENABLE(fortran,
     [  --enable-fortran    Turn on fortran code snipets],
     [case "${enableval}" in
       yes) fortran=true ;;
       no)  fortran=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-fortran) ;;
     esac],[fortran=false])

AC_ARG_ENABLE(emass,
     [  --enable-emass    Turn on electronic mass optimization],
     [case "${enableval}" in
       yes) emass=true ;;
       no)  emass=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-emass) ;;
     esac],[emass=false])

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

AC_ARG_ENABLE(molecularity,
     [  --enable-molecularity    Turn on molecularity],
     [case "${enableval}" in
       yes) molecularity=true ;;
       no)  molecularity=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-molecularity) ;;
     esac],[molecularity=true])

AC_ARG_ENABLE(galib,
     [  --enable-galib    Turn on galib],
     [case "${enableval}" in
       yes) galib=true ;;
       no)  galib=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-galib) ;;
     esac],[galib=false])

AC_ARG_ENABLE(ce,
     [  --enable-cedarwin    Turn on ce darwin],
     [case "${enableval}" in
       yes) ce=true ;;
       no)  ce=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-ce) ;;
     esac],[ce=true])

AC_ARG_ENABLE(pdarwin,
     [  --enable-pdarwin    Turn on pescan-darwin],
     [case "${enableval}" in
       yes) pdarwin=true ;;
       no)  pdarwin=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-pdarwin) ;;
     esac],[pdarwin=true])

AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=true])
AH_TEMPLATE([_PRINT_ALL_PROCS])

AH_TEMPLATE([_DOFORTRAN_])
if test x$fortran = xtrue || test x$emass = xtrue; then
  AC_PROG_FC
  NREL_AC_PROG_FCPP
  AC_FC_LIBRARY_LDFLAGS
  # turn on automatic name mangling
  AC_FC_WRAPPERS
fi
if test x$fortran = xtrue; then
  AC_DEFINE([_DOFORTRAN_],[],[])
fi

AX_COMPILER_VENDOR

# enables vff and pescan for GA
if test x$pdarwin = xtrue || test x$molecularity = xtrue || test $emass = true; then
  pescan=true
  vff=true
fi
if test x$pdarwin = xtrue || test x$molecularity = xtrue|| test x$ce = xtrue || test x$emass = true;
then
  galib=true
fi
if test x$emass = xtrue; then
  lapack=true
fi

if test x$debug = xtrue; then
  CXXFLAGS=`echo $CXXFLAGS | sed "s/-O2//"`
  CFLAGS=`echo $CFLAGS | sed "s/-O2//"`
  AC_DEFINE([_PRINT_ALL_PROCS],[],[])
  AX_CHECK_COMPILER_FLAGS([-fbounds-check], 
                          [CPP_DEBUG_FLAGS+=" -fbounds-check"],)
else
  CXXFLAGS=`echo $CXXFLAGS | sed "s/-g//"`
  CFLAGS=`echo $CFLAGS | sed "s/-g//"`
  CPP_DEBUG_FLAGS=
fi
AC_SUBST(CPP_DEBUG_FLAGS)
ismultiproc=`grep processor /proc/cpuinfo | wc -l`
if test $ismultiproc > 1; then
  ismultiproc="true"
else
  ismultiproc="false"
fi
AC_ARG_ENABLE(mpi,
[  --enable-mpi    Turn on mpi],
[case "${enableval}" in
  yes) mpi=true ;;
  no)  mpi=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
esac],[mpi=$ismultiproc])
AH_TEMPLATE([_MPI])

AC_ARG_ENABLE(nolaunch,
[  --enable-nolaunch    Runs fake pescan/darwin],
[case "${enableval}" in
  yes) nolaunch=true ;;
  no)  nolaunch=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-nolaunch) ;;
esac],[nolaunch=false])
AH_TEMPLATE([_NOLAUNCH])
if test x$nolaunch = xtrue; then
  AC_DEFINE([_NOLAUNCH],[],[])
fi

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

# conditional
AM_CONDITIONAL(DOMPI, test x$mpi = xtrue)
AM_CONDITIONAL(DOPESCAN, test x$pescan = xtrue)
AM_CONDITIONAL(DOMOLECULARITY, test x$molecularity = xtrue)
AM_CONDITIONAL(DOCE, test x$ce = xtrue)
AM_CONDITIONAL(DOVFF, test x$vff = xtrue)
AM_CONDITIONAL(DOPDARWIN, test x$pdarwin = xtrue)
AM_CONDITIONAL(NOLAUNCH, test x$nolaunch = xtrue)
AM_CONDITIONAL(DODARWIN, test x$galib = xtrue)
AM_CONDITIONAL(DOLAPACK, test x$lapack = xtrue)
AM_CONDITIONAL(DOEMASS, test x$emass = xtrue)
AM_CONDITIONAL(DOFORTRAN, test x$fortran = xtrue)

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

#  eo libs iff galib is enabled 
if test x$galib = xtrue; then
  PKG_CHECK_MODULES(EO, eo >= 1.0.1 )
fi 
AC_SUBST(GSL_INCLUDE)
AC_SUBST(GSL_LIB)

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

# checks for STLs
AC_CXX_HAVE_STL
if test "$ac_cv_cxx_have_stl" = yes; then
  LIBS=`echo ${LIBS} | sed "s/-lstdc++//" `
  LIBS="${LIBS} -lm "
  case "$ax_cv_c_compiler_vendor" in 
    gnu) LIBS=" ${LIBS} -lstdc++";;
    intel) LIBS=" ${LIBS} -lstdc++ -lcxa -lcxaguard";;
    ibm) FFLAGS="-qsuffix=f=f90 -lC";;
    portland) LIBS=" ${LIBS} -L/usr/lib -lstdc++ ";;
  esac
 else
   AC_MSG_ERROR("Could not find standard template libraries")
fi


AC_OUTPUT(Makefile atat/Makefile physics/Makefile opt/Makefile lamarck/Makefile
          vff/Makefile mpi/Makefile pescan_interface/Makefile darwin/Makefile 
          print/Makefile lapack/Makefile )

echo ""; echo ""
echo "Will Compile following modules: "
echo " _ CE related stuff in lamarck directory"

if test $galib = true; then
  echo " _ GA libraries"
fi
if test $vff = true; then
  echo " _ vff functional ";
fi
if test $pescan = true; then
  echo " _ interface to nanopse/pescan functional ";
fi
if test $pdarwin = true; then
  echo " _ Genetic Algorithm for Band Gap"
fi
if test $ce = true; then
  echo " _ Genetic Algorithm for Cluster Expansion"
fi
if test x$molecularity = xtrue; then
  echo " _ Genetic Algorithm for Molecularity"
fi
if test x$emass = xtrue; then
  echo " _ Genetic Algorithm for Electronic Mass Optimization"
fi
echo ""
if test $mpi = true; then
  echo "Will Include MPI where available"
fi
if test $debug = true; then
  echo "Compiling with debug options"
fi
if test $nolaunch = true; then
  echo "Compiling with _NOLAUNCH.... Are you sure you know what you're doing? "
fi
echo ""; echo ""
