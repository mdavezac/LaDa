#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @NUMPY_CXXFLAGS@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(math_python_h) 
pkgpyexec_LTLIBRARIES = math.la 
math_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
math_la_SOURCES = $(math_python_cc) 
math_la_LDFLAGS = $(AM_LDFLAGS) -module
math_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
