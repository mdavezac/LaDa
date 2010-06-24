#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
              @BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @EIGEN_CXXFLAGS@
AM_LDFLAGS =

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(vff_python_h) 
subpkgpyexec_PYTHON = __init__.py
subpkgpyexec_LTLIBRARIES = vff.la
subpkgpyexecdir = $(pkgpyexecdir)/vff
vff_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
vff_la_SOURCES = $(vff_python_cc) 
vff_la_LDFLAGS = $(AM_LDFLAGS) -module
vff_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
