bin_PROGRAMS = vff
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
vff_SOURCES = main.cc 
vff_CPPFLAGS =  $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
vff_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@ @FCLIBS@
vff_LDADD = ../lamarck/liblamarck.la ../atat/libatat.la libvff.la\
            ../physics/libphysics.la @MPI_ADDLIB@ 
if DOFORTRAN 
vff_LDADD += ../opt/libopt.la 
endif
if DOFORTRAN
  noinst_PROGRAMS = fortran_test
  fortran_test_LD = $(FC) 
  fortran_test_SOURCES = fortran_test.f90 
  fortran_test_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
  fortran_test_LDADD = lib_vff.la 
  lib_LTLIBRARIES = lib_vff.la
  lib_vff_la_SOURCES = fortran.cc functional.cc 
  lib_vff_la_LIBADD = ../physics/libphysics.la ../lamarck/liblamarck.la ../atat/libatat.la \
  		    ../print/libprint.la @MPI_ADDLIB@
  lib_vff_ladir = $(includedir)/vff
  nobase_lib_vff_la_HEADERS = fortran.h
endif
noinst_LTLIBRARIES = libvff.la
libvff_la_SOURCES = functional.cc layered.cc va.cc
libvff_ladir = $(includedir)/vff
nobase_libvff_la_HEADERS = functional.h layered.h va.h
libvff_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
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
