#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(MINIMIZER_PYTHON_H) 
lib_LTLIBRARIES = libpyminimizer.la 
libpyminimizer_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpyminimizer_la_SOURCES = $(MINIMIZER_PYTHON_CC) 
libpyminimizer_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libpyminimizer_la_LIBADD =  ../libminimizer.la 
../libminimizer.la:
	cd ../; $(MAKE) libminimizer.la
endif
