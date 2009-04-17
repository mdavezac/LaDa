#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_PYTHON@ @BOOST_SERIALIZATION@ 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(CRYSTAL_PYTHON_H) 
pkgpyexec_LTLIBRARIES = crystal.la 
crystal_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
crystal_la_SOURCES = $(CRYSTAL_PYTHON_CC) 
crystal_la_LDFLAGS = $(AM_LDFLAGS) -module
crystal_la_LIBADD =  ../../python/libpython_lada.la 
../../python/libpython_lada.la:
	cd ../../python; $(MAKE) libpython_lada.la
endif
