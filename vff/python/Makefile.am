#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(VFF_PYTHON_H) 
pkgpyexec_LTLIBRARIES = vff.la 
vff_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
vff_la_SOURCES = $(VFF_PYTHON_CC) 
vff_la_LDFLAGS = $(AM_LDFLAGS) -module
vff_la_LIBADD =  ../../python/libpython_lada.la 
../../python/libpython_lada.la:
	cd ../../python; $(MAKE) libpython_lada.la
endif
