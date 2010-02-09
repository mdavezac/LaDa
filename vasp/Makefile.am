if DOPYTHON
include $(top_srcdir)/Files

SUBDIRS = extract incar methods

subpkgpyexec_LTLIBRARIES = _vasp.la 
subpkgpyexec_PYTHON = __init__.py mesh.py launch.py kpoints.py specie.py files.py specie.py
subpkgpyexecdir = $(pkgpyexecdir)/vasp

_vasp_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
_vasp_la_SOURCES = vasp.cc
_vasp_la_LDFLAGS = $(AM_LDFLAGS) -module
_vasp_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
endif
