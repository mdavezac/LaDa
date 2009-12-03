#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(pescan_python_h) 
subpkgpyexec_PYTHON = __init__.py va.py
subpkgpyexecdir = $(pkgpyexecdir)/escan
subpkgpyexec_LTLIBRARIES = _escan.la 
_escan_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
_escan_la_SOURCES = $(pescan_python_cc) 
_escan_la_LDFLAGS = $(AM_LDFLAGS) -module 
_escan_la_LIBADD =  $(top_builddir)/python/libpython_lada.la
endif
