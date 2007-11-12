bin_PROGRAMS = vff fortran_test
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
vff_SOURCES = main.cc 
vff_CPPFLAGS =  $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
vff_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
vff_LDADD = ../lamarck/liblamarck.la ../atat/libatat.la libvff.la\
	             ../physics/libphysics.la @MPI_ADDLIB@
fortran_test_SOURCES = fortran_test.f90 
fortran_test_CPPFLAGS =  $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
fortran_test_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
fortran_test_LDADD = lib_vff.la @MPI_ADDLIB@
lib_LTLIBRARIES = lib_vff.la
lib_vff_la_SOURCES = fortran.cc functional.cc 
lib_vff_la_LIBADD = ../physics/libphysics.la ../lamarck/liblamarck.la ../atat/libatat.la \
	            ../print/libprint.la @MPI_ADDLIB@ 
noinst_LTLIBRARIES = libvff.la
libvff_la_SOURCES = functional.cc layered.cc
libvff_ladir = vff
nobase_libvff_la_HEADERS = functional.h layered.h
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
