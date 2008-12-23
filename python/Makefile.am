#
# Version: $Id$
#
include $(top_srcdir)/Files
if DOCE
  CECC=$(PYTHON_CE_CC)
  CEH=$(PYTHON_CE_H)
  CELIB=../ce/libce.la
endif 
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@\
	      @MINUIT_CFLAGS@ @GSL_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@  @BOOST_LDFLAGS@ @MINUIT_LIBS@ @GSL_LIBS@ @BOOST_PYTHON@ @BOOST_FS@  
lib_LTLIBRARIES = liblada.la 
liblada_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
liblada_la_SOURCES = $(PYTHON_CC) $(CECC)
liblada_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
liblada_la_LIBADD =  ../crystal/libcrystal.la \
		     ../atat/libatat.la \
		     ../opt/libpythonconvenience.la \
		     ../minimizer/libminimizer.la \
		     ../physics/libphysics.la \
		     ../vff/libvff.la \
                     @MPI_ADDLIB@  
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

../vff/libvff.la:
	cd ../vff; $(MAKE) libvff.la

../minimizer/libminimizer.la:
	cd ../minimizer; $(MAKE) libminimizer.la

../opt/libpythonconvenience.la : 
	cd ../opt; $(MAKE) libpythonconvenience.la

../opt/errors.o : ../opt/errors.cc
	cd ../opt; $(MAKE) errors.o
