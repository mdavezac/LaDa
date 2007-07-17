bin_PROGRAMS = vff
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
vff_SOURCES = main.cc 
vff_CPPFLAGS =  $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
vff_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
vff_LDADD = ../lamarck/liblamarck.la ../atat/libatat.la libvff.la ../physics/libphysics.la @MPI_ADDLIB@
noinst_LTLIBRARIES = libvff.la
libvff_la_SOURCES = functional.cc
libvff_ladir = vff
nobase_libvff_la_HEADERS = functional.h
libvff_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
	cd ../lamarck; $(MAKE)
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
