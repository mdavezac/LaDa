#
# Version: $Id$
#
include ../Files
if DOCE
  CECC=$(PYTHON_CE_CC)
  CEH=$(PYTHON_CE_H)
  CELIB=../ce/libce.la
endif 
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @TINYXML_LIBS@  @BOOST_LDFLAGS@ @BOOST_PYTHON@
lib_LTLIBRARIES = liblada.la 
liblada_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) @BOOST_CXXFLAGS@ \
              	-I@PYTHON_INCLUDE_DIR@
liblada_la_SOURCES = $(PYTHON_CC) $(CECC)
liblada_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
liblada_la_LIBADD =  $(CELIB) \
                     ../crystal/libcrystal.la \
		     ../atat/libatat.la \
		     ../physics/libphysics.la \
                     @MPI_ADDLIB@ 
../crystal/licrystal.la: 
	cd ../crystal; $(MAKE) 

../atat/libatat.la: 
	cd ../atat; $(MAKE) 

../atat/libphysic.la: 
	cd ../physics; $(MAKE) 

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)

../ce/libce.la:
	cd ../ce; $(MAKE)
