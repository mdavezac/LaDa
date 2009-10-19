#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @NUMPY_CXXFLAGS@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(ce_python_h) 
pkgpyexec_LTLIBRARIES = ce.la 
ce_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) -D _CUBIC_CE_
ce_la_SOURCES = $(ce_python_cc) 
ce_la_LDFLAGS = $(AM_LDFLAGS) -module
ce_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
