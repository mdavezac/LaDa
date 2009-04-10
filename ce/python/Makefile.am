#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(CE_PYTHON_H) 
pkgpyexec_LTLIBRARIES = ce.la 
ce_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) -D _CUBIC_CE_
ce_la_SOURCES = $(CE_PYTHON_CC) 
ce_la_LDFLAGS = $(AM_LDFLAGS) -module
models_la_LIBADD =  ../../python/libpython_lada.la 
../../python/libpython_lada.la:
	cd ../../python; $(MAKE) libpython_lada.la
endif
