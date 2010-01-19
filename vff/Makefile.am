#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPI_CFLAGS@ @GSL_CFLAGS@ @BOOST_CXXFLAGS@ @MINUIT_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @GSL_LIBS@ @BOOST_FS@  @BOOST_SYSTEM@ \
	     @BOOST_MPI@ @BOOST_SERIALIZATION@ @FCLIBS@ @MINUIT_LIBS@
fortranSOURCE = 
fortranHEADER = 
# fortranSOURCE += $(VFF_FORTRAN_CC)
# fortranHEADER += $(VFF_FORTRAN_H)

bin_PROGRAMS = vff layered_vff
# noinst_PROGRAMS = $(fortranBIN)

noinst_LTLIBRARIES = libvff.la
libvff_la_SOURCES = $(vff_cc) $(fortranSOURCE)
libvff_ladir = $(includedir)/LaDa/vff
nobase_libvff_la_HEADERS = $(vff_h) $(fortranHEADER)
libvff_la_CPPFLAGS = $(AM_CPPFLAGS) 

vff_SOURCES = main.cc 
vff_CPPFLAGS =  $(AM_CPPFLAGS) 
vff_LDFLAGS =  $(AM_LDFLAGS) @BOOST_PROGOP@ 
vff_LDADD = libvff.la \
            $(top_builddir)/crystal/libcrystal.la \
	    $(top_builddir)/math/libmath.la \
            $(top_builddir)/physics/libphysics.la \
            $(top_builddir)/print/libprint.la \
            $(top_builddir)/opt/libopt.la \
	    $(top_builddir)/minimizer/libminimizer.la \
	    $(fortranLDADD) \
	    @MPI_ADDLIB@ \
            $(fortranLDFLAGS) \
            @BOOST_MPI@ @BOOST_SERIALIZATION@
layered_vff_SOURCES = $(vff_SOURCES)
layered_vff_CPPFLAGS =  $(vff_CPPFLAGS) -D _LAYERED
layered_vff_LDFLAGS =  $(vff_LDFLAGS)
layered_vff_LDADD = $(vff_LDADD)

SUBDIRS = python
# fortran_test_SOURCES = fortran_test.f90 
# fortran_test_LDFLAGS =  $(vff_LDFLAGS)
# fortran_test_LDADD = $(vff_LDADD)
