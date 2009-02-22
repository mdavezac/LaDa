#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PYTHON@ @MINUIT_LIBS@ @GSL_LIBS@ @FCLIBS@ \
	     @BOOST_FS@ @BOOST_SERIALIZATION@ @FCLIBS@ @LD_IAGA_FLAGS@ 

EXTRA_DIST = $(PESCAN_PYTHON_H) 
pkgpyexec_LTLIBRARIES = escan.la 
escan_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
escan_la_SOURCES = $(PESCAN_PYTHON_CC) 
escan_la_LDFLAGS = $(AM_LDFLAGS) -module 
escan_la_LIBADD =  ../libescan.la \
                   ../../atat/libatat.la \
                   ../../print/libprint.la \
                   ../../minimizer/libminimizer.la \
                   ../../crystal/libcrystal.la \
                   ../../opt/libpythonconvenience.la \
                   ../../vff/libvff.la
../libescan.la:
	cd ../; $(MAKE) libescan.la
../../crystal/libcrystal.la:
	cd ../; $(MAKE) libcrystal.la
../../atat/libatat.la:
	cd ../../atat; $(MAKE) libatat.la
../../opt/libpythonconvenience.la:
	cd ../../opt; $(MAKE) libpythonconvenience.la
../../print/libprint.la:
	cd ../../print; $(MAKE) libprint.la
../../vff/libvff.la:
	cd ../../vff; $(MAKE) libvff.la
../../minimizer/libminimizer.la:
	cd ../../minimizer; $(MAKE) libminimizer.la
endif
