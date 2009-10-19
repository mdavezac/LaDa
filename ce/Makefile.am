#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @BOOST_CXXFLAGS@ @MPI_CFLAGS@ -I $(top_srcdir) 
AM_LDFLAGS = @TINYXML_LIBS@
if DOCE
bin_PROGRAMS = regular cubic tetra 
cubic_SOURCES = main.cc 
cubic_CPPFLAGS =  $(AM_CPPFLAGS) -D_CUBIC_CE_
cubic_LDFLAGS =  $(AM_LDFLAGS)  @BOOST_LDFLAGS@ \
                 @BOOST_MPI@ @BOOST_SERIALIZATION@ @BOOST_PROGOP@ @BOOST_FS@ \
		 @BOOST_SYSTEM@ @BOOST_REGEX@ 
cubic_LDADD =   libce.la \
                $(top_builddir)/crystal/libcrystal.la \
		$(top_builddir)/atat/libatat.la \
		$(top_builddir)/physics/libphysics.la \
		@MPI_ADDLIB@ \
		$(top_builddir)/print/libprint.la \
                $(top_builddir)/opt/libopt.la \
		$(top_builddir)/minimizer/libminimizer.la
tetra_SOURCES = $(cubic_SOURCES)
tetra_CPPFLAGS =  $(AM_CPPFLAGS) -D_TETRAGONAL_CE_
tetra_LDFLAGS = $(cubic_LDFLAGS)
tetra_LDADD = $(cubic_LDADD) 
regular_SOURCES = rmain.cc
regular_CPPFLAGS =  $(AM_CPPFLAGS)  -D_CUBIC_CE_ @GSL_CFLAGS@ @MINUIT_CFLAGS@
regular_LDFLAGS = $(cubic_LDFLAGS) @GSL_LIBS@ @MINUIT_LIBS@
regular_LDADD = libce.la \
                $(top_builddir)/crystal/libcrystal.la \
		$(top_builddir)/atat/libatat.la \
		$(top_builddir)/physics/libphysics.la \
		$(top_builddir)/print/libprint.la \
                $(top_builddir)/opt/libopt.la \
		$(top_builddir)/minimizer/libminimizer.la \
		$(top_builddir)/mpi/libmpiobject.la
endif
noinst_LTLIBRARIES = libce.la 
libce_la_SOURCES = $(ce_cc)
libce_ladir = $(includedir)/LaDa/ce
nobase_libce_la_HEADERS = $(ce_h) $(ce_impl_h)
libce_la_CPPFLAGS =  $(AM_CPPFLAGS) 

SUBDIRS = python
