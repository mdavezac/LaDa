#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @MINUIT_CFLAGS@ @GSL_CFLAGS@
AM_LDFLAGS = 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(MINIMIZER_PYTHON_H) 
pkgpyexec_LTLIBRARIES = minimizer.la 
minimizer_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
minimizer_la_SOURCES = $(MINIMIZER_PYTHON_CC) 
minimizer_la_LDFLAGS = $(AM_LDFLAGS) -module
minimizer_la_LIBADD =  ../../python/libpython_lada.la 
../../python/libpython_lada.la:
	cd ../../python; $(MAKE) libpython_lada.la
endif
