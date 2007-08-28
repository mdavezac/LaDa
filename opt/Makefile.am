AM_CPPFLAGS = @MPI_DEFS@ @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPICH_CFLAGS@ @GSL_INCLUDE@ 
AM_LDFLAGS = @TINYXML_LIBS@
# noinst_LTLIBRARIES = libopt.la
libopt_ladir = opt
# libopt_la_SOURCES = dbrent.F90  df1dim.F90  frprmn.F90  linmin.F90
nobase_libopt_la_HEADERS = bisection.h ndim_iterator.h opt_minimize_base.h opt_minimize_opt++.h \
			   opt_polynome.h compose_functors.h  opt_frprmn.h  opt_minimize_gsl.h \
		       	   opt_minimize_ssquared.h  types.h convex_hull.h opt_function_base.h \
			   opt_minimize.h opt_minimize_wang.h va_minimizer.h fitness_function.h \
			   opt_functors.h opt_minimize_linear.h opt_monome.h traits.h
libopt_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 

