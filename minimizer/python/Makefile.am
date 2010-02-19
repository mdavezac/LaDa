if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @MINUIT_CFLAGS@ @GSL_CFLAGS@\
       	@NUMPY_CXXFLAGS@ @EIGEN_CXXFLAGS@
AM_LDFLAGS = 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(minimizer_python_h) 
pkgpyexec_LTLIBRARIES = minimizer.la 
minimizer_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
minimizer_la_SOURCES = $(minimizer_python_cc) 
minimizer_la_LDFLAGS = $(AM_LDFLAGS) -module
minimizer_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
