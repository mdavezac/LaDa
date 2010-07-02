#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPI_CFLAGS@ @GSL_CFLAGS@ @BOOST_CXXFLAGS@ @MINUIT_CFLAGS@\
              @EIGEN_CXXFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @GSL_LIBS@ @BOOST_FS@  @BOOST_SYSTEM@ \
	     @BOOST_MPI@ @BOOST_SERIALIZATION@ @FCLIBS@ @MINUIT_LIBS@

bin_PROGRAMS = vff layered_vff

noinst_LTLIBRARIES = libvff.la
libvff_la_SOURCES = $(vff_cc) 
libvff_ladir = $(includedir)/LaDa/vff
nobase_libvff_la_HEADERS = $(vff_h) 
libvff_la_CPPFLAGS = $(AM_CPPFLAGS) 

vff_SOURCES = main.cc 
vff_CPPFLAGS =  $(AM_CPPFLAGS)  @GSL_CFLAGS@
vff_LDFLAGS =  $(AM_LDFLAGS) @BOOST_PROGOP@ 
vff_LDADD = libvff.la \
            $(top_builddir)/crystal/libcrystal.la \
	    $(top_builddir)/math/libmath.la \
            $(top_builddir)/physics/libphysics.la \
            $(top_builddir)/opt/libopt.la \
	    $(top_builddir)/minimizer/libminimizer.la \
	    @MPI_ADDLIB@ \
            @BOOST_MPI@ @BOOST_SERIALIZATION@
layered_vff_SOURCES = $(vff_SOURCES)
layered_vff_CPPFLAGS =  $(vff_CPPFLAGS) -D _LAYERED
layered_vff_LDFLAGS =  $(vff_LDFLAGS)
layered_vff_LDADD = $(vff_LDADD)

SUBDIRS = python
