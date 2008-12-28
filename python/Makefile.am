#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @TINYXML_LIBS@  @BOOST_LDFLAGS@ @BOOST_PYTHON@ @BOOST_FS@  \
             @BOOST_MPI@ @BOOST_SERIALIZATION@ @BOOST_MPI_PYTHON@

if DOCE
  CECC=$(PYTHON_CE_CC)
  CEH=$(PYTHON_CE_H)
  CELIB=../ce/libce.la
endif 

if DOPYTHONPESCAN
  ESCANCC=$(PYTHON_ESCAN_CC)
  ESCANH=$(PYTHON_ESCAN_H)
  ESCANLIB= ../mpi/libmpiobject.la ../pescan_interface/libescan.la \
	    ../print/libprint.la @LD_IAGA_FLAGS@
  AM_CPPFLAGS += @MINUIT_CFLAGS@ @GSL_CFLAGS@
  AM_LDFLAGS += @FCLIBS@ 
endif

if DOVFF
  VFFCC=$(PYTHON_VFF_CC)
  VFFH=$(PYTHON_VFF_H)
  VFFLIB=../minimizer/libminimizer.la ../vff/libvff.la
  AM_CPPFLAGS += @MINUIT_CFLAGS@ @GSL_CFLAGS@
  AM_LDFLAGS += @MINUIT_LIBS@ @GSL_LIBS@ 
endif

lib_LTLIBRARIES = liblada.la 
liblada_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
liblada_la_SOURCES = $(PYTHON_CC) $(CECC) $(ESCANCC) $(VFFCC) 
liblada_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
liblada_la_LIBADD =  $(CELIB) $(VFFLIB) $(ESCANLIB)\
                     ../crystal/libcrystal.la \
		     ../atat/libatat.la \
		     ../opt/libpythonconvenience.la \
		     ../physics/libphysics.la \
                     @MPI_ADDLIB@  
if DOPYTHONPESCAN
  liblada_la_FCFLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  liblada_la_F77FLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  liblada_la_F90FLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  liblada_la_FFLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
endif
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
	
../pescan_interface/libescan.la: 
	cd ../pescan_interface; $(MAKE) libescan.la

../print/libprint.la: 
	cd ../print; $(MAKE) libprint.la
