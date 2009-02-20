#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@ 

EXTRA_DIST = $(OPT_PYTHON_H) 
lib_LTLIBRARIES = libpyopt.la 
libpyopt_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpyopt_la_SOURCES = $(OPT_PYTHON_CC) 
libpyopt_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libpyopt_la_LIBADD =  ../libpythonconvenience.la 
../libpythonconvenience.la:
	cd ../; $(MAKE) libpythonconvenience.la
EXTRA_DIST = $(OPT_PYTHON_H) 
pkgpyexec_LTLIBRARIES = _opt.la 
_opt_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
_opt_la_SOURCES = $(OPT_PYTHON_CC) 
_opt_la_LDFLAGS = $(AM_LDFLAGS) -module
_opt_la_LIBADD =  ../libpythonconvenience.la 
../libpythonconvenience.la:
	cd ../; $(MAKE) libpythonconvenience.la
endif
