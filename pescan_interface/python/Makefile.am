#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(PESCAN_PYTHON_H) 
lib_LTLIBRARIES = libpyescan.la 
libpyescan_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpyescan_la_SOURCES = $(PESCAN_PYTHON_CC) 
libpyescan_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libpyescan_la_LIBADD =  ../libescan.la 
../libescan.la:
	cd ../; $(MAKE) libescan.la
endif
