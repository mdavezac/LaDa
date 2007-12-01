bin_PROGRAMS = escan
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ 
if NOLAUNCH
  AM_CPPFLAGS += @EO_CFLAGS@
  AM_LDFLAGS  += @EO_LIBS@
endif    
escan_SOURCES = main.cc
escan_CPPFLAGS = $(AM_CPPFLAGS)  -I ..
escan_LDFLAGS =  $(AM_LDFLAGS) @MPICH_LIBS@
escan_LDADD = libpescan.la ../lamarck/liblamarck.la ../atat/libatat.la ../vff/libvff.la\
              ../physics/libphysics.la ../print/libprint.la @MPI_ADDLIB@ 
if DOFORTRAN
  escan_LDFLAGS +=  @FCLIBS@
  escan_LDADD += ../opt/libopt.la 
else
  escan_CPPFLAGS += @GSL_INCLUDE@
  escan_LDFLAGS += @GSL_LIB@
endif
noinst_LTLIBRARIES = libpescan.la
libpescan_la_SOURCES = interface.cc va.cc
libpescan_ladir = $(includedir)/pescan_interface 
nobase_libpescan_la_HEADERS = interface.h va.h
libpescan_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
	cd ../lamarck; $(MAKE)
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
../vff/libvff.la:
	cd ../vff; $(MAKE)
../print/libprint.la:
	cd ../print; $(MAKE)
if DOFORTRAN
../mpi/libopt.la:
	cd ../opt; $(MAKE)
endif
