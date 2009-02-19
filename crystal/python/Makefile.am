#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(CRYSTAL_PYTHON_H) 
lib_LTLIBRARIES = libpycrystal.la 
libpycrystal_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpycrystal_la_SOURCES = $(CRYSTAL_PYTHON_CC) 
libpycrystal_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libpycrystal_la_LIBADD =  ../libcrystal.la 
../libcrystal.la:
	cd ../; $(MAKE) libcrystal.la
endif
