include ../Files
AM_CPPFLAGS = @MPI_DEFS@ @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPICH_CFLAGS@ @GSL_INCLUDE@ 
AM_LDFLAGS = @TINYXML_LIBS@
if DOFORTRAN
noinst_LTLIBRARIES = libopt.la
libopt_la_SOURCES = $(OPT_F90)
libopt_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
endif

libopt_ladir = $(includedir)/opt

nobase_libopt_la_HEADERS = $(OPT_H)
