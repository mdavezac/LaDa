#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PYTHON@ @MINUIT_LIBS@ @GSL_LIBS@ @FCLIBS@ \
	     @BOOST_FS@ @BOOST_SERIALIZATION@

EXTRA_DIST = $(VFF_PYTHON_H) 
pkgpyexec_LTLIBRARIES = vff.la 
vff_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
vff_la_SOURCES = $(VFF_PYTHON_CC) 
vff_la_LDFLAGS = $(AM_LDFLAGS) -module
vff_la_LIBADD =  ../libvff.la \
		 ../../atat/libatat.la \
		 ../../print/libprint.la \
		 ../../minimizer/libminimizer.la \
		 ../../crystal/libcrystal.la \
		 ../../opt/libpythonconvenience.la
../libvff.la:
	cd ../; $(MAKE) libvff.la
../../crystal/libcrystal.la:
	cd ../; $(MAKE) libcrystal.la
../../atat/libatat.la:
	cd ../../atat; $(MAKE) libatat.la
../../opt/libpythonconvenience.la:
	cd ../../opt; $(MAKE) libpythonconvenience.la
../../minimizer/libminimizer.la:
	cd ../../minimizer; $(MAKE) libminimizer.la
../../print/libprint.la:
	cd ../../print; $(MAKE) libprint.la
endif
