#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @MINUIT_CFLAGS@ @GSL_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PYTHON@ @MINUIT_LIBS@ @GSL_LIBS@ @FCLIBS@ \
	     @BOOST_FS@
if DOMPI
  AM_LDFLAGS += @BOOST_SERIALIZATION@ @BOOST_MPI@
endif

EXTRA_DIST = $(MINIMIZER_PYTHON_H) 
pkgpyexec_LTLIBRARIES = minimizer.la 
minimizer_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
minimizer_la_SOURCES = $(MINIMIZER_PYTHON_CC) 
minimizer_la_LDFLAGS = $(AM_LDFLAGS) -module
minimizer_la_LIBADD =  ../libminimizer.la \
		       ../../opt/libpythonconvenience.la
../libminimizer.la:
	cd ../; $(MAKE) libminimizer.la
../../opt/libpythonconvenience.la:
	cd ../../opt; $(MAKE) libpythonconvenience.la
endif
