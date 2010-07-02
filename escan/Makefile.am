if DOPYTHON 
AM_CPPFLAGS = -I$(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
              @BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@ @NUMPY_CXXFLAGS@\
              @EIGEN_CXXFLAGS@ 
AM_LDFLAGS = @BOOST_LDFLAGS@ @BOOST_PYTHON@ 

install-exec-local: ../../python/libpython_lada.la

EXTRA_DIST = escan.hpp bandgap.hpp emass.hpp dipole_elements.hpp wavefunctions.hpp
subpkgpyexec_PYTHON = __init__.py va.py derivatives.py band_structure.py _extract.py \
		      _bandgap.py _wfns.py _methods.py
subpkgpyexecdir = $(pkgpyexecdir)/escan
subpkgpyexec_LTLIBRARIES = _escan.la 
_escan_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
_escan_la_SOURCES = escan.cc python.cc wavefunctions.cc escan_types.f90 \
		    wavefunctions_module.f90 fortran.f90
_escan_la_LDFLAGS = $(AM_LDFLAGS) -module 
_escan_la_FCFLAGS = @IAGA_FCFLAGS@
endif 
