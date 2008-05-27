AC_INIT([LaDa],[3.2])
AM_INIT_AUTOMAKE()
# fi

AC_PROG_CC
AC_PROG_CXX

AC_PROG_EGREP
dompidefault=false
if [[[ $(echo $CXX | $EGREP "mpi(cxx|c\+\+|CC)" - ) == $CXX ]]]; then 
 dompidefault=true
 mpicompilers=true
fi

AC_PROG_INSTALL

AM_CONFIG_HEADER(config.h)

AC_ARG_ENABLE(all,
     [  --enable-all                  Turn off/on all (most) modules],
     [case "${enableval}" in
       yes) enall=true ;;
       no)  enall=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-all) ;;
     esac],[enall=true])

AC_ARG_ENABLE(force-mpi,
     [  --enable-force-mpi             Forces mpi to be on],
     [case "${enableval}" in
       yes) force_mpi=true ;;
       no)  force_mpi=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-force-mpi) ;;
     esac],[force_mpi=$dompidefault])
AC_ARG_ENABLE(mpi,
     [  --enable-mpi                  Turn on mpi],
     [case "${enableval}" in
       yes) mpi=true ;;
       no)  mpi=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
     esac],[mpi=$force_mpi])
AH_TEMPLATE([_MPI])

AC_ARG_ENABLE(fortran,
     [  --enable-fortran              Turn on fortran minimization.
                                Incompatible with gsl/nag.],
     [case "${enableval}" in
       yes) fortran=true ;;
       no)  fortran=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-fortran) ;;
     esac],[fortran=false])

AC_ARG_ENABLE(nag,
     [  --enable-nag                  Turn on NAG minimization.
                                Incompatible with gsl/fortran],
     [case "${enableval}" in
       yes) nag=true ;;
       no)  nag=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-fortran) ;;
     esac],[nag=false])

AC_ARG_ENABLE(emass,
     [  --enable-emass                Turn on electronic mass optimization],
     [case "${enableval}" in
       yes) emass=true ;;
       no)  emass=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-emass) ;;
     esac],[emass=$enall])

AC_ARG_ENABLE(vff,
     [  --enable-vff                  Turn on vff],
     [case "${enableval}" in
       yes) vff=true ;;
       no)  vff=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-vff) ;;
     esac],[vff=$enall])

AC_ARG_ENABLE(pescan,
     [  --enable-pescan               Turn on pescan],
     [case "${enableval}" in
       yes) pescan=true ;;
       no)  pescan=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-pescan) ;;
     esac],[pescan=$enall])

AC_ARG_ENABLE(molecularity,
     [  --enable-molecularity         Turn on molecularity],
     [case "${enableval}" in
       yes) molecularity=true ;;
       no)  molecularity=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-molecularity) ;;
     esac],[molecularity=$enall])

AC_ARG_ENABLE(galib,
     [  --enable-galib                Turn on galib],
     [case "${enableval}" in
       yes) galib=true ;;
       no)  galib=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-galib) ;;
     esac],[galib=$enall])

AC_ARG_ENABLE(ceopt,
     [  --enable-ceopt                Turn on cluster-expansion optimization],
     [case "${enableval}" in
       yes) ceopt=true ;;
       no)  ceopt=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-ceopt) ;;
     esac],[ceopt=$enall])

AC_ARG_ENABLE(bandgapopt,
     [  --enable-bandgapopt           Turn on bandgap optimization],
     [case "${enableval}" in
       yes) bandgapopt=true ;;
       no)  bandgapopt=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-bandgap) ;;
     esac],[bandgapopt=$enall])

AC_ARG_ENABLE(emassopt,
     [  --enable-emassopt             Turn on electronic effective mass optimization],
     [case "${enableval}" in
       yes) emassopt=true ;;
       no)  emassopt=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-emassopt) ;;
     esac],[emassopt=$enall])

AC_ARG_ENABLE(debug,
     [  --enable-debug                Turn on debugging],
     [case "${enableval}" in
       yes) debug=true ;;
       no)  debug=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
     esac],[debug=true])
AH_TEMPLATE([_PRINT_ALL_PROCS])

AC_ARG_ENABLE(nolaunch,
     [  --enable-nolaunch             Runs fake pescan/darwin],
     [case "${enableval}" in
       yes) nolaunch=true ;;
       no)  nolaunch=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-nolaunch) ;;
     esac],[nolaunch=false])
AH_TEMPLATE([_NOLAUNCH])
if test x$nolaunch = xtrue; then
  AC_DEFINE([_NOLAUNCH],[],[])
fi

AC_ARG_ENABLE(iaga-dir,
              [  --enable-iaga-dir=?  specifies the base nanopse
                                directory, eg $HOME/usr/src/nanopse ],
              [ iaga_dir=${enableval} ] )
AC_ARG_ENABLE(iaga-ld,
              [  --enable-iaga-ld=?  specifies the base iaga ld flags ],
              [ iaga_ld=${enableval} ] )

AC_ARG_ENABLE(boost-ld,
              [  --enable-boost-ld=?  specifies the base boost ld flags ],
              [ boost_ld=${enableval} ] )
AC_ARG_ENABLE(boost-cxx,
              [  --enable-boost-cxx=?  specifies the boost cxx flags ],
              [ boost_cxx=${enableval} ] )


if test x$bandgapopt = xtrue ; then 
  galib=true
  pescan=true
  vff=true
fi
if test x$molecularity = xtrue; then
  galib=true
  pescan=true
  vff=true
fi
if test x$emassopt = xtrue; then
  galib=true
  emass=true
  pescan=true
  vff=true
fi
if test x$ceopt = xtrue; then
  galib=true
fi
lapack=false
if test x$emass = xtrue; then
  pescan=true
  vff=true
  lapack=true
fi

AH_TEMPLATE([_DONAG])
AH_TEMPLATE([_DOGSL])
if test $vff = false; then
  gsl=false
  fortran=false
  nag=false
fi
if test $nag = true; then
  gsl=false
  fortran=false
  AC_DEFINE([_DONAG],[],[])
  PKG_CHECK_MODULES(NAG, nag )
fi
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
if test $fortran = true; then
  gsl=false
  nag=false
else
  gsl=true
fi
AC_SUBST(NAG_CFLAGS)
AC_SUBST(NAG_LIBS)
if test x$gsl = xtrue; then
  fortran=false
  AC_DEFINE([_DOGSL],[],[])
  PKG_CHECK_MODULES(GSL, gsl >= 1.8 )
fi 

# ************************************
# ****** Checks compiler compliances
# ************************************

AH_TEMPLATE([_DOFORTRAN])
fortrancompiler=false
if test x$fortran = xtrue; then
 fortrancompiler=true
fi
if test x$lapack = xtrue; then
  fortrancompiler=true
fi
if test $pescan = true; then
  fortrancompiler=true;
fi
if test x$fortrancompiler = xtrue; then
  AC_PROG_FC
  NREL_AC_PROG_FCPP
  AC_FC_LIBRARY_LDFLAGS
  # turn on automatic name mangling
  AC_FC_WRAPPERS
fi
if test x$fortran = xtrue; then
  AC_DEFINE([_DOFORTRAN],[],[])
fi

AH_TEMPLATE([_DIRECTIAGA])
direct_iaga=false;
if test "x$pescan$mpi" = "xtruetrue"; then
  if test "x$iaga_dir" = "x" ; then
    direct_iaga=false;
  else
    direct_iaga=true;
#   FC_IAGA_FLAGS="$FCFLAGS -J$iaga_dir/pescan/lib -J$iaga_dir/genpot/lib"
     
    cp -u $iaga_dir/genpot/lib/mpigroup.mod pescan_interface
    cp -u $iaga_dir/genpot/lib/MPIGROUP.mod pescan_interface
    cp -u $iaga_dir/pescan/lib/escan_comp_api.mod pescan_interface
    cp -u $iaga_dir/pescan/lib/ESCAN_COMP_API.mod pescan_interface
    cp -u $iaga_dir/pescan/lib/eigenenergy.mod pescan_interface
    cp -u $iaga_dir/pescan/lib/EIGENENERGY.mod pescan_interface
    IAGA_DIR="$iaga_dir/"
    LD_IAGA_FLAGS=" $iaga_dir/pescan/lib/libpescan.a  $iaga_dir/genpot/lib/libgenpot.a "
    AC_DEFINE([_DIRECTIAGA],[1],[])
    if test "x$iaga_ld" = "x" ; then
      LD_IAGA_FLAGS=" $iaga_dir/pescan/lib/libpescan.a  $iaga_dir/genpot/lib/libgenpot.a "
    else
      LD_IAGA_FLAGS="$LD_IAGA_FLAGS $iaga_ld"
    fi
  fi
fi
AC_SUBST(LD_IAGA_FLAGS )
AC_SUBST(IAGA_DIR)
# AC_SUBST(FC_IAGA_FLAGS )

AX_COMPILER_VENDOR

# C++ compiler compliance
AC_CXX_FULL_SPECIALIZATION_SYNTAX
AC_CXX_PARTIAL_SPECIALIZATION
AC_CXX_DEFAULT_TEMPLATE_PARAMETERS
AC_CXX_FULL_DEFAULT_TEMPLATE_PARAMETERS

# autoscan suggestions
AC_CXX_HAVE_IEEE_MATH
AC_CXX_COMPLEX_MATH_IN_NAMESPACE_STD
AC_CXX_HAVE_COMPLEX_MATH1
AC_CHECK_FUNCS([gettimeofday],[],[{(exit 1); exit 1; }] )
AC_CHECK_HEADERS([float.h],[],[{(exit 1); exit 1; }] )
AC_CHECK_HEADERS([limits.h],[],[{(exit 1); exit 1; }] )
AC_CHECK_HEADERS([sys/time.h],[],[{(exit 1); exit 1; }] )
AC_HEADER_STDBOOL
AC_HEADER_TIME

AH_TEMPLATE([_DEBUG])
if test x$debug = xtrue; then
  CXXFLAGS=`echo $CXXFLAGS | sed "s/-O2//"`
  CFLAGS=`echo $CFLAGS | sed "s/-O2//"`
  AC_DEFINE([_PRINT_ALL_PROCS],[],[])
# AX_CHECK_COMPILER_FLAGS([-fbounds-check], 
#                         [CPP_DEBUG_FLAGS+=" -fbounds-check"],)
  AC_DEFINE([_DEBUG],[1],[])
else
  CXXFLAGS=`echo $CXXFLAGS | sed "s/-g//"`
  CFLAGS=`echo $CFLAGS | sed "s/-g//"`
  CPP_DEBUG_FLAGS=
fi
AC_SUBST(CPP_DEBUG_FLAGS)

# ************************************
# ****** End of Checks compiler compliances
# ************************************


#  MPI libs
if test x$mpi = xtrue; then
  if test $force_mpi = true; then
    AC_DEFINE([_MPI],[1],[])
    MPI_ADDLIB="../mpi/libmpiobject.la"
    MPI_CFLAGS="$MPI_CFLAGS -D _OPENMPI_MPI_"
    mpicompilers=true

  else
    mpicompilers=false
    PKG_CHECK_MODULES(WHATEVER, mpich >= 1.2.7, [mpichisthere=yes], [mpichisthere=no] )
    PKG_CHECK_MODULES(WHATEVER, openmpi >= 1.2, [openmpiisther=yes], [openmpiisther=no] )
    if test $openmpiisther = yes; then
      PKG_CHECK_MODULES(MPI, openmpi >= 1.2 )
      AC_DEFINE([_MPI],[1],[])
      MPI_ADDLIB="../mpi/libmpiobject.la"
      MPI_CFLAGS="$MPI_CFLAGS -D _OPENMPI_MPI_"
    elif test $mpichisthere = yes; then
      PKG_CHECK_MODULES(MPI, mpich >= 1.2.7 )
      AC_DEFINE([_MPI],[1],[])
      MPI_ADDLIB="../mpi/libmpiobject.la"
      MPI_CFLAGS="$MPI_CFLAGS -D _MPICH_MPI_"
    fi
  fi
fi
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_LIBS)
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
AM_CONDITIONAL(DOLAYEREDOPT, test x$molecularity = xtrue)
AM_CONDITIONAL(DOCEOPT, test x$ceopt = xtrue)
AM_CONDITIONAL(DOVFF, test x$vff = xtrue)
AM_CONDITIONAL(DOBANDGAPOPT, test x$bandgapopt = xtrue)
AM_CONDITIONAL(NOLAUNCH, test x$nolaunch = xtrue)
AM_CONDITIONAL(DODARWIN, test x$galib = xtrue)
AM_CONDITIONAL(DOLAPACK, test x$lapack = xtrue)
AM_CONDITIONAL(DOEMASS, test x$emass = xtrue)
AM_CONDITIONAL(DOEMASSOPT, test x$emassopt = xtrue)
AM_CONDITIONAL(DOFORTRAN, test x$fortran = xtrue)
AM_CONDITIONAL(DOGSL, test x$gsl = xtrue)
AM_CONDITIONAL(DONAG, test x$nag = xtrue)
AM_CONDITIONAL(DODIRECTIAGA, test x$direct_iaga = xtrue)


#  eo libs iff galib is enabled 
if test x$galib = xtrue; then
  PKG_CHECK_MODULES(EO, eo >= 1.0.1 )
fi 
AC_SUBST(NAG_CFLAGS)
AC_SUBST(NAG_LIBS)

# checks for compiler specific flags 
#
COMPILER_CXXFLAGS=
# first checks if friend injection is needed (for atat) (gcc >= 4.0)
# if test "$ac_compiler_gnu"="yes"; then
# AX_CHECK_COMPILER_FLAGS([-malign-double], 
#                         [COMPILER_CXXFLAGS+=" -malign-double"],)
# fi
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
    portland) LIBS=" ${LIBS} -lC ";;
  esac
 else
   AC_MSG_ERROR("Could not find standard template libraries")
fi

AC_CHECK_PROG(founddox, doxygen, true, false)
AM_CONDITIONAL(DODOC, test x$founddox = xtrue)
AC_CHECK_PROG(FOUNDDOT, dot, "YES", "NO")

DOCDIR="doc"
if test $datadir = "\${prefix}/share"; then
  DOCDIR="doc"
else
  DOCDIR="$datadir"
fi
AC_SUBST(DOCDIR)

# boost stuff
if test "x$boost_ld" = "x" ; then
  BOOST_LDFLAGS="-L/usr/lib -lboost_signals -lboost_program_options"
else
  BOOST_LDFLAGS="$boost_ld"
fi
if test "x$boost_cxx" = "x" ; then
  BOOST_CXXFLAGS="-I /usr/include"
else
  BOOST_CXXFLAGS="$boost_cxx"
fi
AC_SUBST(BOOST_LDFLAGS)
AC_SUBST(BOOST_CXXFLAGS)

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
if test $pescan$direct_iaga = truefalse; then
  echo " _ interface to nanopse/pescan functional as external program. ";
fi
if test $pescan$direct_iaga = truetrue; then
  echo " _ interface to nanopse/pescan functional as external library. ";
  echo "   iaga base directory = $iaga_dir ";
fi
if test x$emass = xtrue; then
  echo " _ Effective Electronic Mass calculations for ternary (100) uperlattices" 
fi
if test $bandgapopt = true; then
  echo " _ Genetic Algorithm for Band Gap"
fi
if test $ceopt = true; then
  echo " _ Genetic Algorithm for Cluster Expansion"
fi
if test x$molecularity = xtrue; then
  echo " _ Genetic Algorithm for Molecularity"
fi
if test x$emassopt = xtrue; then
  echo " _ Genetic Algorithm for Electronic Mass Optimization"
fi
echo ""
if test $fortran = true; then
  echo "Will use fortran minimizer"
fi
if test $mpi$mpichisthere = trueyes; then
  echo "Will Include MPI where available, using mpich"
elif test $mpi$openmpiisther = trueyes; then
  echo "Will Include MPI where available, using openmpi"
elif test $mpi$mpicompilers = truetrue; then
  echo "Will Include MPI where available, using openmpi"
fi
if test $lapack = true; then
  echo "Will Include some lapack/blas routines for Matrix diagonalization"
fi
if test $debug = true; then
  echo "Compiling with debug options"
fi
if test $nolaunch = true; then
  echo "Compiling with _NOLAUNCH.... Are you sure you know what you're doing? "
fi
echo ""
if test $founddox = true; then
  echo "Type \"make doc\" for doxygen documentation"
  echo "Then check out $DOCDIR/html/index.html"
else
  echo "You will need to install doxygen to generate the documentation"
fi
if test $FOUNDDOT = NO; then
  echo "Install dot for Documentation pretty graph with GraphView"
fi
echo ""; echo ""
