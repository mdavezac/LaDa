#
# Version: $Id$
#
include ../Files
# if DOPYTHON
#   pythonlib = libcrystalpy.la
# endif
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@
AM_LDFLAGS =
noinst_LTLIBRARIES = libcrystal.la # $(pythonlib)
libcrystal_la_SOURCES = $(CRYSTAL_CC)
libcrystal_ladir = $(includedir)/crystal
nobase_libcrystal_la_HEADERS = $(CRYSTAL_H) 
libcrystal_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
# if DOPYTHON
# libcrystalpy_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) @BOOST_CXXFLAGS@ \
#       		     -I@PYTHON_INCLUDE_DIR@
# libcrystalpy_la_SOURCES = $(CE_PYTHON_CC)
# libcrystalpy_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
# libcrystalpy_la_LIBADD = @BOOST_LDFLAGS@ @BOOST_PYTHON@
# endif
../atat/libatat.la: 
	cd ../atat; $(MAKE) 

../atat/libphysic.la: 
	cd ../physics; $(MAKE) 

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)

../print/libprint.la:
	cd ../print; $(MAKE)
