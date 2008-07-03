#
# Version: $Id$
#
include ../Files
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ @BOOST_CXXFLAGS@ \
	      @TINYXML_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
if DOFORTRAN
OPTF90 = $(OPT_F90)
endif
noinst_LTLIBRARIES = libopt.la
libopt_la_SOURCES = $(OPTF90) $(OPT_CC)
libopt_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS)
libopt_ladir = $(includedir)/opt
nobase_libopt_la_HEADERS = $(OPT_H)
