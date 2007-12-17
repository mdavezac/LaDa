include ../Files
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
if DOFORTRAN 
  fortranLDADD   = ../opt/libopt.la 
  fortranLDFLAGS = @FCLIBS@
  fortranBIN     = fortran_test
  fortranSOURCE = $(VFF_FORTRAN_CC)
  fortranHEADER = $(VFF_FORTRAN_H)
  gslCPPFLAGS =
  gslLDFLAGS  =
else 
  fortranLDADD   = 
  fortranLDFLAGS = 
  fortranBIN = 
  fortranSOURCE = 
  fortranHEADER = 
  gslCPPFLAGS =  @GSL_INCLUDE@ 
  gslLDFLAGS  = @GSL_LIB@
endif

bin_PROGRAMS = vff layered_vff
noinst_PROGRAMS = $(fortranBIN)

noinst_LTLIBRARIES = libvff.la
libvff_la_SOURCES = $(VFF_CC) $(fortranSOURCE)
libvff_ladir = $(includedir)/vff
nobase_libvff_la_HEADERS = $(VFF_H) $(fortranHEADER)
libvff_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 

vff_SOURCES = main.cc 
vff_CPPFLAGS =  $(AM_CPPFLAGS) -I .. $(gslCPPFLAGS)
vff_LDFLAGS =  $(AM_LDFLAGS) @MPICH_LIBS@ $(gslLDFLAGS) $(fortranLDFLAGS)
vff_LDADD = ../lamarck/liblamarck.la \
	    ../atat/libatat.la \
	    libvff.la \
            ../physics/libphysics.la \
	    $(fortran_LDADD) \
	    @MPI_ADDLIB@ 
layered_vff_SOURCES = $(vff_SOURCES)
layered_vff_CPPFLAGS =  $(vff_CPPFLAGS) -D _LAYERED
layered_vff_LDFLAGS =  $(vff_LDFLAGS)
layered_vff_LDADD = $(vff_LDADD)

fortran_test_LD = $(FC) 
fortran_test_SOURCES = fortran_test.f90 
fortran_test_LDFLAGS =  $(vff_LDFLAGS)
fortran_test_LDADD = $(vff_LDADD)

../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
	cd ../lamarck; $(MAKE)
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
../mpi/libprint.la:
	cd ../print; $(MAKE)

if DOFORTRAN
../mpi/libopt.la:
	cd ../opt; $(MAKE)
endif
