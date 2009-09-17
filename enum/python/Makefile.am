#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_PYTHON@ 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(enum_python_h) enum.py
pkgpyexec_LTLIBRARIES = enumeration.la 
enumeration_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
enumeration_la_SOURCES = $(enum_python_cc) 
enumeration_la_LDFLAGS = $(AM_LDFLAGS) -module
enumeration_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
