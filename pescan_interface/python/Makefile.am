#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(PESCAN_PYTHON_H) 
pkgpyexec_LTLIBRARIES = escan.la 
escan_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
escan_la_SOURCES = $(PESCAN_PYTHON_CC) 
escan_la_LDFLAGS = $(AM_LDFLAGS) -module 
escan_la_LIBADD =  ../../python/libpython_lada.la
../../python/libpython_lada.la:
	cd ../../python/libpython_lada.la; $(MAKE) libpython_lada.la
endif
