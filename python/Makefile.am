#
# Version: $Id$
#
include ../Files
if DOCE
  CECC=$(PYTHON_CE_CC)
  CEH=$(PYTHON_CE_H)
  CELIB=../ce/libce.la
endif 
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @TINYXML_LIBS@  @BOOST_LDFLAGS@
lib_LTLIBRARIES = liblada.la 
liblada_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
liblada_la_SOURCES = $(PYTHON_CC) $(CECC)
liblada_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
liblada_la_LIBADD =  $(CELIB) \
                     ../crystal/libcrystal.la \
		     ../atat/libatat.la \
		     ../physics/libphysics.la \
		     ../opt/libpythonconvenience.la \
                     @MPI_ADDLIB@  @BOOST_PYTHON@ @BOOST_FS@ 
../crystal/libcrystal.la: 
	cd ../crystal; $(MAKE) libcrystal.la

../atat/libatat.la: 
	cd ../atat; $(MAKE) libatat.la

../atat/libphysic.la: 
	cd ../physics; $(MAKE) libphysics.la

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE) libmpiobject.la

../ce/libce.la:
	cd ../ce; $(MAKE) libce.la

../opt/libpythonconvenience.la : 
	cd ../opt; $(MAKE) libpythonconvenience.la

../opt/errors.o : ../opt/errors.cc
	cd ../opt; $(MAKE) errors.o
