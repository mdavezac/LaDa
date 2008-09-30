#
# Version: $Id$
#
include ../Files
# if DOPYTHON
#   pythonlib = libcrystalpy.la
# endif
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@
AM_LDFLAGS = @BOOST_REGEX@
noinst_LTLIBRARIES = libcrystal.la # $(pythonlib)
libcrystal_la_SOURCES = $(CRYSTAL_CC)
libcrystal_ladir = $(includedir)/LaDa/crystal
nobase_libcrystal_la_HEADERS = $(CRYSTAL_H) 
libcrystal_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 

../atat/libphysic.la: 
	cd ../physics; $(MAKE) 

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)

../print/libprint.la:
	cd ../print; $(MAKE)
