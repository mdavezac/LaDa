#
# Version: $Id$
#
include ../Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPI_CFLAGS@ @GSL_CFLAGS@ @BOOST_CXXFLAGS@ @MINUIT_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @GSL_LIBS@ @MINUIT_LIBS@ \
	     @BOOST_PROGOP@ @BOOST_REGEX@ @BOOST_FS@ \
	     @BOOST_SYSTEM@ @BOOST_MPI@ @BOOST_SERIALIZATION@
if NOLAUNCH
  AM_CPPFLAGS += @EO_CFLAGS@
  AM_LDFLAGS  += @EO_LIBS@
endif    

LDADD = $(top_builddir)/opt/libopt.la 
if DODIRECTIAGA
   LDADD +=  $(top_builddir)/mpi/libmpiobject.la
   AM_LDFLAGS += @FCLIBS@
endif

bin_PROGRAMS = escan layered_escan

noinst_LTLIBRARIES = libescan.la
libescan_ladir = $(includedir)/LaDa/pescan_interface
libescan_la_SOURCES = $(pescan_cc) 
nobase_libescan_la_HEADERS = $(pescan_h) $(pescan_impl_h) 
libescan_la_CPPFLAGS = $(AM_CPPFLAGS) 
if DODIRECTIAGA
  libescan_la_SOURCES += $(pescan_f90)
  libescan_la_FCFLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  libescan_la_F77FLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  libescan_la_F90FLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  libescan_la_FFLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
endif

escan_SOURCES = main.cc
escan_CPPFLAGS = $(AM_CPPFLAGS) 
escan_LDFLAGS =  $(AM_LDFLAGS) 
escan_LDADD = libescan.la \
	      $(top_builddir)/vff/libvff.la \
	      $(top_builddir)/crystal/libcrystal.la \
	      $(top_builddir)/atat/libatat.la \
	      $(top_builddir)/minimizer/libminimizer.la \
              $(top_builddir)/physics/libphysics.la \
	      $(top_builddir)/print/libprint.la \
	      @MPI_ADDLIB@ \
	      @LD_IAGA_FLAGS@ $(LDADD) 


layered_escan_SOURCES = main.cc
layered_escan_CPPFLAGS = $(escan_CPPFLAGS) -D _LAYERED_
layered_escan_LDFLAGS =  $(escan_LDFLAGS) 
layered_escan_LDADD = $(escan_LDADD)

SUBDIRS = python

if DODIRECTIAGA
$(top_builddir)/lapack/liblapack.la: \
    $(addprefix $(top_srcdir)/lapack/, $(lapack_h) )\
    $(addprefix $(top_srcdir)/lapack/, $(lapack_f90) )\
    $(addprefix $(top_srcdir)/lapack/, $(lapack_cc) ) 
	@cd $(top_builddir)/lapack; $(MAKE) liblapack.la
endif






