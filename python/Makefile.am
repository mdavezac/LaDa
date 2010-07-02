#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @BOOST_LDFLAGS@ @BOOST_PYTHON@ @TINYXML_LIBS@ @BOOST_FS@ @BOOST_SERIALIZATION@ 
if DOMPI
  AM_LDFLAGS += @BOOST_MPI@ 
endif

EXTRA_DIST = $(python_h)
lib_LTLIBRARIES = libpython_lada.la 
libpython_lada_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
libpython_lada_la_SOURCES = 
libpython_lada_la_LDFLAGS = $(AM_LDFLAGS)
libpython_lada_la_LIBADD = 
pythonextensions = opt physics math crystal ga
if DOCE
  pythonextensions += ce
  libpython_lada_la_LIBADD += $(top_builddir)/ce/libce.la 
endif 
if DOVFF
  pythonextensions += vff
  libpython_lada_la_LIBADD += $(top_builddir)/vff/libvff.la 
endif
if DOPOTENTIALS
  pythonextensions += potentials
  libpython_lada_la_LIBADD += $(top_builddir)/atomic_potentials/libpotentials.la 
endif

if DOMODELS
  pythonextensions += models 
  libpython_lada_la_LIBADD += $(top_builddir)/models/libmodels.la 
endif
if DOENUM
  pythonextensions += enumeration
  libpython_lada_la_LIBADD += $(top_builddir)/enum/libenum.la 
endif

if DOMINIMIZER
  pythonextensions += minimizer 
  libpython_lada_la_LIBADD += $(top_builddir)/minimizer/libminimizer.la
  libpython_lada_la_LDFLAGS += @GSL_LIBS@  @MINUIT_LIBS@
  libpython_lada_la_CPPFLAGS += @GSL_CFLAGS@ @MINUIT_CFLAGS@ @FCLIBS@
endif
if DOVASP
  pythonextensions += vasp
endif
libpython_lada_la_LIBADD +=  $(top_builddir)/crystal/libcrystal.la \
                             $(top_builddir)/math/libmath.la \
                             $(top_builddir)/opt/libpythonconvenience.la \
                             $(top_builddir)/physics/libphysics.la \
                             $(top_builddir)/print/libprint.la 


pkgpyexec_PYTHON = __init__.py
force: ;
__init__.py: force
	@list='$(pythonextensions)'; \
	 echo "\"\"\" Lamarck-Darwin extension library.\"\"\"" > __init__.py;\
	 for library in $$list; do \
	   echo "import $$library" >> __init__.py; \
	 done
clean-local:
	-rm -f __init__.py

# EXTRA_DIST = $(PYTHON_H) $(PYTHON_CE_H) $(PYTHON_VFF_H) $(PYTHON_ESCAN_H)
$(top_builddir)/opt/errors.o : $(top_srcdir)/opt/errors.cc
	cd $(top_builddir)/opt; $(MAKE) errors.o
