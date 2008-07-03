#
# Version: $Id$
#
include ../Files
if DOPYTHON
  pythonlib = liblamarckpy.la
endif
bin_PROGRAMS = cubic tetra
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
cubic_SOURCES = main.cc 
cubic_CPPFLAGS =  $(AM_CPPFLAGS) -I .. -D_CUBIC_CE_
cubic_LDFLAGS =  $(AM_LDFLAGS)  @MPI_LIBS@ @BOOST_LDFLAGS@ @BOOST_MPI@ @BOOST_PROGOP@
cubic_LDADD = liblamarck.la \
		../atat/libatat.la \
		../physics/libphysics.la \
		@MPI_ADDLIB@ \
		../print/libprint.la 
tetra_SOURCES = $(cubic_SOURCES)
tetra_CPPFLAGS =  $(AM_CPPFLAGS) -I .. -D_TETRAGONAL_CE_
tetra_LDFLAGS = $(cubic_LDFLAGS)
tetra_LDADD = $(cubic_LDADD) 
noinst_LTLIBRARIES = liblamarck.la $(pythonlib)
liblamarck_la_SOURCES = $(CRYSTAL_CC) $(CE_CC)
liblamarck_ladir = $(includedir)/lamarck
nobase_liblamarck_la_HEADERS = $(CRYSTAL_H) $(CE_H)
liblamarck_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
if DOPYTHON
liblamarckpy_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) @BOOST_CXXFLAGS@ \
                           -I@PYTHON_INCLUDE_DIR@
liblamarckpy_la_SOURCES = $(CE_PYTHON_CC)
liblamarckpy_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
liblamarckpy_la_LIBADD = @BOOST_LDFLAGS@ @BOOST_PYTHON@
endif
../atat/libatat.la: 
	cd ../atat; $(MAKE) 

../atat/libphysic.la: 
	cd ../physics; $(MAKE) 

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)

../print/libprint.la:
	cd ../print; $(MAKE)
