#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(CE_PYTHON_H) 
lib_LTLIBRARIES = libpyce.la 
libpyce_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpyce_la_SOURCES = $(CE_PYTHON_CC) 
libpyce_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libpyce_la_LIBADD =  ../libce.la 
../libce.la:
	cd ../; $(MAKE) libce.la
endif
