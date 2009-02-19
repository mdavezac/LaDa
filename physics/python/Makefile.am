#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(PHYSICS_PYTHON_H) 
lib_LTLIBRARIES = libpyphysics.la 
libpyphysics_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpyphysics_la_SOURCES = $(PHYSICS_PYTHON_CC) 
libpyphysics_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libpyphysics_la_LIBADD =  ../libphysics.la 
../libphysics.la:
	cd ../; $(MAKE) libphysics.la
endif
