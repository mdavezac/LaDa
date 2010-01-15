#
#  Version: $Id$
#
AC_INIT([LaDa],[3.30])
AM_INIT_AUTOMAKE([-Wno-portability])

AC_CONFIG_MACRO_DIR([m4])


AC_PROG_CC
AC_PROG_CXX
AM_PATH_PYTHON(2.3)

AC_PROG_EGREP
dompidefault=false
if [[[ $(echo $CXX | $EGREP "mpi(cxx|c\+\+|CC)" - ) == $CXX ]]]; then 
 dompidefault=true
fi

if [[ -e ".env" ]]; then
  . "./.env"
else
  echo "Environment variables in file .env not found."
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

AC_ARG_ENABLE(models,
     [  --enable-models               Compiles Coulom + Lennard-Jhones wrappers.],
     [case "${enableval}" in
       yes) models=true ;;
       no)  models=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-models) ;;
     esac],[models=false])

AC_ARG_ENABLE(potentials,
     [  --enable-potentials           Compiles generalized atomic potentials.],
     [case "${enableval}" in
       yes) potentials=true ;;
       no)  potentials=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-potentials) ;;
     esac],[potentials=false])

AC_ARG_ENABLE(enum,
     [  --enable-enum                 Compiles enumeration library.],
     [case "${enableval}" in
       yes) enum=true ;;
       no)  enum=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-enum) ;;
     esac],[enum=false])

AC_ARG_ENABLE(python,
     [  --enable-python               Compiles python bindings.],
     [case "${enableval}" in
       yes) python=true ;;
       no)  python=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-python) ;;
     esac],[python=false])


AC_ARG_ENABLE(mpi,
     [  --enable-mpi                  Turn on mpi],
     [case "${enableval}" in
       yes) mpi=true ;;
       no)  mpi=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
     esac],[mpi=$dompidefault])
AH_TEMPLATE([_MPI])

# AC_ARG_ENABLE(emass,
#      [  --enable-emass                Turn on electronic mass optimization],
#      [case "${enableval}" in
#        yes) emass=true ;;
#        no)  emass=false ;;
#        *) AC_MSG_ERROR(bad value ${enableval} for --enable-emass) ;;
#      esac],[emass=$enall])

AC_ARG_ENABLE(separable,
     [  --enable-separable            Turn on sum of separable functions module.],
     [case "${enableval}" in
       yes) separable=true ;;
       no)  separable=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-separable) ;;
     esac],[separable=$enall])

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

AC_ARG_ENABLE(galib,
     [  --enable-galib                Turn on galib],
     [case "${enableval}" in
       yes) galib=true ;;
       no)  galib=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-galib) ;;
     esac],[galib=$enall])

AC_ARG_ENABLE(ce,
     [  --enable-ce                   Turn on cluster-expansion module],
     [case "${enableval}" in
       yes) ce=true ;;
       no)  ce=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-ce) ;;
     esac],[ce=$enall])

AC_ARG_ENABLE(ceopt,
     [  --enable-ceopt                Turn on cluster-expansion optimization],
     [case "${enableval}" in
       yes) ceopt=true ;;
       no)  ceopt=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-ceopt) ;;
     esac],[ceopt=$enall])

AC_ARG_ENABLE(alloylayersopt,
     [  --enable-alloylayersopt       Turn on alloylayers configuration space optimization],
     [case "${enableval}" in
       yes) alloylayersopt=true ;;
       no)  alloylayersopt=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-alloylayers) ;;
     esac],[alloylayersopt=$enall])

AC_ARG_ENABLE(debug,
     [  --enable-debug                Turn on debugging],
     [case "${enableval}" in
       yes) debug=true ;;
       no)  debug=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
     esac],[debug=true])
AH_TEMPLATE([_PRINT_ALL_PROCS])
AH_TEMPLATE([NDEBUG])

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

# if test x$emassopt = xtrue; then
#   galib=true
#   emass=true
#   pescan=true
#   vff=true
# fi
if test x$alloylayersopt = xtrue; then
  galib=true
  pescan=true
  vff=true
fi
if test x$ceopt = xtrue; then
  ce=true
  galib=true
fi
direct_iaga=true;

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AH_TEMPLATE([__DOCE])
if test x$ce = xtrue; then
  AC_DEFINE([__DOCE],[],[])
fi
AH_TEMPLATE([__DOVFF])
if test x$vff = xtrue; then
  AC_DEFINE([__DOVFF],[],[])
fi
AH_TEMPLATE([__DOPYTHON])
if test x$python = xtrue; then
  AC_DEFINE([__DOPYTHON],[],[])
fi
AH_TEMPLATE([__DOPYTHONESCAN])
if test x$pescan$python$direct_iaga = xtruetruetrue; then
  AC_DEFINE([__DOPYTHONESCAN],[],[])
fi

# ************************************
# ****** Checks compiler compliances
# ************************************

AC_PROG_F77
AC_PROG_FC
NREL_AC_PROG_FCPP
AC_FC_LIBRARY_LDFLAGS
# turn on automatic name mangling
AC_FC_WRAPPERS

AH_TEMPLATE([_DIRECTIAGA])
AC_DEFINE([_DIRECTIAGA],[1],[])
AC_SUBST(LD_IAGA_FLAGS)
AC_SUBST(IAGA_FCFLAGS)
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

AH_TEMPLATE([_LADADEBUG])
if test x$debug = xtrue; then
  CXXFLAGS=`echo $CXXFLAGS | sed "s/-O2//"`
  CFLAGS=`echo $CFLAGS | sed "s/-O2//"`
  AC_DEFINE([_PRINT_ALL_PROCS],[],[])
# AX_CHECK_COMPILER_FLAGS([-fbounds-check], 
#                         [CPP_DEBUG_FLAGS+=" -fbounds-check"],)
  AC_DEFINE([_LADADEBUG],[1],[])
else
  CXXFLAGS=`echo $CXXFLAGS | sed "s/-g//"`
  CFLAGS=`echo $CFLAGS | sed "s/-g//"`
  CPP_DEBUG_FLAGS=
  AC_DEFINE([NDEBUG],[1],[])
fi
CXXFLAGS="$CXXFLAGS"
AC_SUBST(CPP_DEBUG_FLAGS)

# ************************************
# ****** End of Checks compiler compliances
# ************************************


#  MPI libs
if test x$mpi = xtrue; then
  AC_DEFINE([_MPI],[1],[])
  MPI_ADDLIB="../mpi/libmpiobject.la"
  MPI_CFLAGS="-D _OPENMPI_MPI_"
else
  MPI_CFLAGS=""
  MPI_ADDLIB=""
fi
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_ADDLIB)

#  XML libs
AC_SUBST(TINYXML_LIBS)
AC_SUBST(TINYXML_CFLAGS)

#  libm 
AC_CHECK_LIB(m,main)

# conditional
AM_CONDITIONAL(DOMPI, test x$mpi = xtrue)
AM_CONDITIONAL(DOPESCAN, test x$pescan = xtrue)
AM_CONDITIONAL(DOCE, test x$ce = xtrue)
AM_CONDITIONAL(DOCEOPT, test x$ceopt = xtrue)
AM_CONDITIONAL(DOVFF, test x$vff = xtrue)
AM_CONDITIONAL(NOLAUNCH, test x$nolaunch = xtrue)
AM_CONDITIONAL(DODARWIN, test x$galib = xtrue)
AM_CONDITIONAL(DOALLOYLAYERSOPT, test x$alloylayersopt = xtrue)
AM_CONDITIONAL(DODIRECTIAGA, test x$direct_iaga = xtrue)
AM_CONDITIONAL(DOPYTHONPESCAN, test x$python$pescan$direct_iaga = xtruetruetrue)
AM_CONDITIONAL(DOPYTHON, test x$python = xtrue)
AM_CONDITIONAL(DOSEPARABLE, test x$separable = xtrue)
AM_CONDITIONAL(DOMODELS, test x$models = xtrue)
AM_CONDITIONAL(DOENUM, test x$enum = xtrue)
AM_CONDITIONAL(DOPOTENTIALS, test x$potentials = xtrue)
AM_CONDITIONAL(DOMINIMIZER, test xtrue = xtrue)


#  eo libs iff galib is enabled 
if test x$galib != xtrue; then
  EO_CFLAGS=""
  EO_LIBS=""
fi 
AC_SUBST(EO_CFLAGS)
AC_SUBST(EO_LIBS)

#
AC_PROG_LIBTOOL

# checks for STLs
AC_CXX_HAVE_STL
if test "$ac_cv_cxx_have_stl" = yes; then
  LIBS=`echo ${LIBS} | sed "s/-lstdc++//" `
  LIBS="${LIBS} -lm "
  case "$ax_cv_c_compiler_vendor" in 
    gnu) LIBS=" ${LIBS} -lstdc++";;
#   intel) LIBS=" ${LIBS} -lstdc++ -lcxa -lcxaguard";;
    ibm) FFLAGS="-qsuffix=f=f90 -lC";;
    portland)
     LIBS=" ${LIBS} -lC ";
     AH_TEMPLATE([BOOST_UBLAS_UNSUPPORTED_COMPILER])
     AC_DEFINE([BOOST_UBLAS_UNSUPPORTED_COMPILER],[0],[])
     ;;
  esac
 else
   AC_MSG_ERROR("Could not find standard template libraries")
fi

AC_CHECK_PROG(founddox, doxygen, true, false)
AM_CONDITIONAL(DODOC, test x$founddox = xtrue)
AC_CHECK_PROG(FOUNDDOT, dot, "YES", "NO")

DOCDIR="doc"
if test "$datadir" = "\${prefix}/share"; then
  DOCDIR="doc"
else
  DOCDIR="$datadir"
fi
AC_SUBST(DOCDIR)

# boost stuff
# if test "x$boost_ld" = "x" ; then
#   BOOST_LDFLAGS="-L/usr/lib -lboost_signals -lboost_program_options"
# else
#   BOOST_LDFLAGS="$boost_ld"
# fi
# if test "x$boost_cxx" = "x" ; then
#   BOOST_CXXFLAGS="-I /usr/include"
# else
#   BOOST_CXXFLAGS="$boost_cxx"
# fi
if test "$mpi" != true; then
  BOOST_MPI=""
fi
if test "$python" = true; then
  if test "x$PYTHON_INCLUDE_DIR$PYTHON_LIB" = "x"; then
    AX_PYTHON
  fi
else
  BOOST_PYTHON=""
fi
AC_SUBST(PYTHON_INCLUDE_DIR)
AC_SUBST(PYTHON_LIB)
AC_SUBST(BOOST_LDFLAGS)
AC_SUBST(BOOST_MPI)
AC_SUBST(BOOST_SERIALIZATION)
AC_SUBST(BOOST_PROGOP)
AC_SUBST(BOOST_FS)
AC_SUBST(BOOST_SYSTEM)
AC_SUBST(BOOST_SIGNAL)
AC_SUBST(BOOST_PYTHON)
AC_SUBST(BOOST_MPI_PYTHON)
AC_SUBST(BOOST_REGEX)
AC_SUBST(BOOST_CXXFLAGS)
AC_SUBST(NUMPY_CXXFLAGS)
AC_SUBST(MINUIT_LIBS)
AC_SUBST(MINUIT_CFLAGS)

AC_OUTPUT(Makefile math/Makefile physics/Makefile physics/python/Makefile
          opt/Makefile opt/python/Makefile crystal/Makefile crystal/python/Makefile
          vff/Makefile vff/python/Makefile mpi/Makefile pescan_interface/Makefile darwin/Makefile 
          pescan_interface/python/Makefile
          print/Makefile separable/Makefile ce/Makefile ce/python/Makefile
          python/Makefile separable/fixedlattice/Makefile factory/Makefile
          minimizer/Makefile minimizer/python/Makefile models/Makefile models/python/Makefile 
          atomic_potentials/Makefile atomic_potentials/python/Makefile 
          atomic_potentials/collapse/Makefile enum/Makefile enum/python/Makefile 
          ga/Makefile vasp/Makefile )

echo ""; echo ""
echo "Will Compile following modules: "
echo " _ Crystal module library"

if test "$ce" = true; then
  echo " _ Cluster Expansion library module"
fi
if test "$galib" = true; then
  echo " _ GA libraries"
fi
if test "$vff" = true; then
  echo " _ vff functional ";
fi
if test "$pescan$direct_iaga" = truefalse; then
  echo " _ interface to nanopse/pescan functional as external program. ";
fi
if test "$pescan$direct_iaga" = truetrue; then
  echo " _ interface to nanopse/pescan functional as external library.";
fi
# if test "x$emass" = xtrue; then
#   echo " _ Effective Electronic Mass calculations for ternary (100) uperlattices" 
# fi
if test "$ceopt" = true; then
  echo " _ Genetic Algorithm for Cluster Expansion"
fi
if test "$models" = true; then
  echo " _ Model functionals"
fi
if test "$potentials" = true; then
  echo " _ Generalized atomic potentials."
fi
if test "$enum" = true; then
  echo " _ Enumeration library"
fi
# if test x$emassopt = xtrue; then
#   echo " _ Genetic Algorithm for Electronic Mass Optimization"
# fi
if test x$alloylayersopt = xtrue; then
  echo " _ Genetic Algorithm for Alloy Layers configuration space."
fi
if test x$separable = xtrue; then
  echo " _ Sum of separable functions fitting."
fi
echo ""
if test "$mpi" = true; then
  echo "Will Include MPI where available."
fi
if test "x$python" = "xtrue"; then
  echo "Will compile python bindings."
fi
if test "$debug" = true; then
  echo "Compiling with debug options"
fi
if test "$nolaunch" = true; then
  echo "Compiling with _NOLAUNCH.... Are you sure you know what you're doing? "
fi
echo ""
if test "$founddox" = true; then
  echo "Type \"make doc\" for doxygen documentation"
  echo "Then check out $DOCDIR/html/index.html"
else
  echo "You will need to install doxygen to generate the documentation"
fi
if test "$FOUNDDOT" = NO; then
  echo "Install dot for Documentation pretty graph with GraphView"
fi
echo ""; echo ""
