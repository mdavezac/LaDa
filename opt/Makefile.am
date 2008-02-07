include ../Files
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = 
if DOFORTRAN
noinst_LTLIBRARIES = libopt.la
libopt_la_SOURCES = $(OPT_F90)
libopt_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
endif

libopt_ladir = $(includedir)/opt

nobase_libopt_la_HEADERS = $(OPT_H)
