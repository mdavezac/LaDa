#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@ -I $(top_srcdir) 
AM_LDFLAGS = @TINYXML_LIBS@
if DOCE
bin_PROGRAMS = regular cubic tetra 
cubic_SOURCES = main.cc 
cubic_CPPFLAGS =  $(AM_CPPFLAGS) $(top_srcdir) -D_CUBIC_CE_
cubic_LDFLAGS =  $(AM_LDFLAGS)  @BOOST_LDFLAGS@\
                 @BOOST_MPI@ @BOOST_SERIALIZATION@ @BOOST_PROGOP@ @BOOST_FS@ \
		 @BOOST_SYSTEM@ @BOOST_REGEX@ 
cubic_LDADD =   libce.la \
                ../crystal/libcrystal.la \
		../atat/libatat.la \
		../physics/libphysics.la \
		@MPI_ADDLIB@ \
		../print/libprint.la \
                ../opt/libopt.la \
		../minimizer/libminimizer.la \
		@BOOST_FS@
tetra_SOURCES = $(cubic_SOURCES)
tetra_CPPFLAGS =  $(AM_CPPFLAGS) -D_TETRAGONAL_CE_
tetra_LDFLAGS = $(cubic_LDFLAGS)
tetra_LDADD = $(cubic_LDADD) 
regular_SOURCES = rmain.cc
regular_CPPFLAGS =  $(AM_CPPFLAGS) -D_CUBIC_CE_ @GSL_CFLAGS@ @MINUIT_CFLAGS@
regular_LDFLAGS = $(cubic_LDFLAGS) @GSL_LIBS@ @MINUIT_LIBS@
regular_LDADD = libce.la \
                ../crystal/libcrystal.la \
		../atat/libatat.la \
		../physics/libphysics.la \
		../print/libprint.la \
                ../opt/libopt.la \
		../minimizer/libminimizer.la \
		../mpi/libmpiobject.la
endif
noinst_LTLIBRARIES = libce.la 
libce_la_SOURCES = $(CE_CC)
libce_ladir = $(includedir)/LaDa/ce
nobase_libce_la_HEADERS = $(CE_H) $(CE_IMPL_H)
libce_la_CPPFLAGS =  $(AM_CPPFLAGS) 
../crystal/libcrystal.la: 
	cd ../crystal; $(MAKE) libcrystal.la

../atat/libatat.la: 
	cd ../atat; $(MAKE) libatat.la

../atat/libphysic.la: 
	cd ../physics; $(MAKE) libphysics.la

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE) libmpiobject.la

../print/libprint.la:
	cd ../print; $(MAKE) libprint.la

../opt/libopt.la:
	cd ../opt; $(MAKE) libopt.la

../minimizer/libminimizer.la:
	cd ../minimizer; $(MAKE) libminimizer.la
