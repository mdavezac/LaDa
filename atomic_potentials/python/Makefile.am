#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_PYTHON@ 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(potentials_python_h) allsq.py
pkgpyexec_LTLIBRARIES = potentials.la 
potentials_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
potentials_la_SOURCES = $(potentials_python_cc) 
potentials_la_LDFLAGS = $(AM_LDFLAGS) -module
potentials_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
