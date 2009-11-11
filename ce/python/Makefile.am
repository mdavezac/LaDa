#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @NUMPY_CXXFLAGS@
AM_LDFLAGS =  @BOOST_LDFLAGS@ @BOOST_PYTHON@

EXTRA_DIST = $(ce_python_h) 
subpkgpyexec_PYTHON = __init__.py _fit.py
subpkgpyexec_LTLIBRARIES = _ce.la 
subpkgpyexecdir = $(pkgpyexecdir)/ce
_ce_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) -D _CUBIC_CE_
_ce_la_SOURCES = $(ce_python_cc) 
_ce_la_LDFLAGS = $(AM_LDFLAGS) -module
_ce_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
force: ;
__init__.py: force
	@echo "from _ce import *"  > __init__.py; \
	 echo "from _fit import *"  >> __init__.py; 
clean-local:
	-rm -f __init__.py
endif
