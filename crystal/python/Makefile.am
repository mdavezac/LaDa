#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PYTHON@ @BOOST_SERIALIZATION@ 

EXTRA_DIST = $(CRYSTAL_PYTHON_H) 
pkgpyexec_LTLIBRARIES = crystal.la 
crystal_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
crystal_la_SOURCES = $(CRYSTAL_PYTHON_CC) 
crystal_la_LDFLAGS = $(AM_LDFLAGS) -module
crystal_la_LIBADD =  ../libcrystal.la \
		     ../../atat/libatat.la \
		     ../../print/libprint.la \
		     ../../opt/libpythonconvenience.la
../libcrystal.la:
	cd ../; $(MAKE) libcrystal.la
../../atat/libatat.la:
	cd ../../atat; $(MAKE) libatat.la
../../opt/libpythonconvenience.la:
	cd ../../opt; $(MAKE) libpythonconvenience.la
../../print/libprint.la:
	cd ../../print; $(MAKE) libprint.la
endif
