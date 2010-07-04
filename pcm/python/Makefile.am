#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@


EXTRA_DIST = $(models_python_h) 
pkgpyexec_LTLIBRARIES = models.la 
models_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
models_la_SOURCES = $(models_python_cc) 
models_la_LDFLAGS = $(AM_LDFLAGS) -module
models_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
