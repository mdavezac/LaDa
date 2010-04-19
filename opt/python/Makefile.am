if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @EIGEN_CXXFLAGS@
AM_LDFLAGS =

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = $(opt_python_h) 
subpkgpyexec_PYTHON = __init__.py secant.py tempdir.py changedir.py decorators.py
subpkgpyexec_LTLIBRARIES = _opt.la 
subpkgpyexecdir = $(pkgpyexecdir)/opt
_opt_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
_opt_la_SOURCES = $(opt_python_cc) ../fortran_redirect.f90
_opt_la_LDFLAGS = $(AM_LDFLAGS) -module
_opt_la_LIBADD =  $(top_builddir)/python/libpython_lada.la 
__init__.py: force
	@echo "\"\"\"Miscellaneous\"\"\""  > __init__.py; \
	 echo "from _opt import *"  >> __init__.py; \
	 echo "from _opt import __load_vasp_in_global_namespace__ "  >> __init__.py; 
	 echo "from _opt import __load_pescan_in_global_namespace__ "  >> __init__.py; 
force: ;
endif
