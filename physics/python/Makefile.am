#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(PHYSICS_PYTHON_H) 
pkgpyexec_LTLIBRARIES = physics.la 
physics_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
physics_la_SOURCES = $(PHYSICS_PYTHON_CC) 
physics_la_LDFLAGS = $(AM_LDFLAGS) -module
physics_la_LIBADD =  ../../python/libpython_lada.la 
../../python/libpython_lada.la:
	cd ../../python; $(MAKE) libpython_lada.la
endif
