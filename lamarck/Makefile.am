bin_PROGRAMS = lamarck
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
lamarck_SOURCES = main.cc 
lamarck_CPPFLAGS =  $(AM_CPPFLAGS) -I .. 
lamarck_LDFLAGS =  $(AM_LDFLAGS)  -lstdc++ @MPICH_LIBS@
lamarck_LDADD = liblamarck.la ../atat/libatat.la ../physics/libphysics.la @MPI_ADDLIB@
noinst_LTLIBRARIES = liblamarck.la
liblamarck_la_SOURCES = cluster.cc  constituent_strain.cc functional_builder.cc \
			harmonic.cc  lattice.cc polynome.cc structure.cc 
liblamarck_ladir = lamarck
nobase_liblamarck_la_HEADERS = atom.h cluster.h  lattice.h polynome.h \
                               constituent_strain.h  functional_builder.h  harmonic.h \
			       structure.h
liblamarck_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 

../atat/libphysic.la: 
	cd ../physics; $(MAKE) 

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
