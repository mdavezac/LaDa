AM_CPPFLAGS = @CPP_DEBUG_FLAGS@ @COMPILER_CXXFLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = 
noinst_LTLIBRARIES = libmpiobject.la
libmpiobject_la_SOURCES = mpi_object.cc
libmpiobject_la_HEADERS = mpi_object.h
libmpiobject_ladir = $(includedir)/mpi
libmpiobject_la_LDFLAGS = $(AM_LDFLAGS) @MPICH_LIBS@
libmpiobject_la_CPPFLAGS = $(AM_CPPFLAGS) $(MPICH_CFLAGS)



bin_PROGRAMS = mpitest
mpitest_SOURCES = main.cc 
mpitest_LDFLAGS =  $(AM_LDFLAGS) @MPICH_LIBS@
mpitest_LDADD = @MPI_ADDLIB@
