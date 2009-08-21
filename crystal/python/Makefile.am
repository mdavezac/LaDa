#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_PYTHON@ @BOOST_SERIALIZATION@ 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(crystal_python_h) 
pkgpyexec_LTLIBRARIES = crystal.la 
crystal_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
crystal_la_SOURCES = $(crystal_python_cc) 
crystal_la_LDFLAGS = $(AM_LDFLAGS) -module
crystal_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
