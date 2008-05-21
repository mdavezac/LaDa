include ../Files
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ @BOOST_CXXFLAGS@
noinst_LTLIBRARIES = libprint.la
libprint_ladir = $(includedir)/print
libprint_la_SOURCES = $(PRINT_CC)
nobase_libprint_la_HEADERS = $(PRINT_H)
libprint_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
