#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(CE_PYTHON_H) 
pkgpyexec_LTLIBRARIES = _ce.la 
_ce_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
_ce_la_SOURCES = $(CE_PYTHON_CC) 
_ce_la_LDFLAGS = $(AM_LDFLAGS) -module
_ce_la_LIBADD =  ../libce.la 
../libce.la:
	cd ../; $(MAKE) libce.la
endif
