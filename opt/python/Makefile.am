#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(opt_python_h) 
pkgpyexec_LTLIBRARIES = opt.la 
opt_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
opt_la_SOURCES = $(opt_python_cc) 
opt_la_LDFLAGS = $(AM_LDFLAGS) -module
opt_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
