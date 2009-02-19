#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(ATAT_PYTHON_H) 
lib_LTLIBRARIES = libpyatat.la 
libpyatat_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpyatat_la_SOURCES = $(ATAT_PYTHON_CC) 
libpyatat_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libpyatat_la_LIBADD =  ../libatat.la 
../libatat.la:
	cd ../; $(MAKE) libatat.la
endif
