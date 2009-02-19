#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(VFF_PYTHON_H) 
lib_LTLIBRARIES = libpyvff.la 
libpyvff_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpyvff_la_SOURCES = $(VFF_PYTHON_CC) 
libpyvff_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libpyvff_la_LIBADD =  ../libvff.la 
../libvff.la:
	cd ../; $(MAKE) libvff.la
endif
