if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir)@CPP_DEBUG_FLAGS@ @BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ 
AM_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_PYTHON@ 

SUBDIRS = extract incar methods

subpkgpyexec_LTLIBRARIES = _vasp.la 
subpkgpyexec_PYTHON = __init__.py mesh.py launch.py kpoints.py specie.py files.py specie.py
subpkgpyexecdir = $(pkgpyexecdir)/vasp

_vasp_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
_vasp_la_SOURCES = vasp.cc
_vasp_la_LDFLAGS = $(AM_LDFLAGS) $(VASP_LDFLAGS) -module
_vasp_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
