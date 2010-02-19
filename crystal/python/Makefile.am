#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ @EIGEN_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @NUMPY_CXXFLAGS@
AM_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_PYTHON@ @BOOST_SERIALIZATION@ 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(crystal_python_h) 
subpkgpyexec_PYTHON = __init__.py cellsym.py gruber.py
subpkgpyexec_LTLIBRARIES = _crystal.la 
subpkgpyexecdir = $(pkgpyexecdir)/crystal
_crystal_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
_crystal_la_SOURCES = $(crystal_python_cc) 
_crystal_la_LDFLAGS = $(AM_LDFLAGS) -module
_crystal_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
