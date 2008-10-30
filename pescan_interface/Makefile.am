#
# Version: $Id$
#
include ../Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPI_CFLAGS@ @NAG_CFLAGS@ @GSL_CFLAGS@ @BOOST_CXXFLAGS@ 
AM_LDFLAGS = @TINYXML_LIBS@ @NAG_LIBS@ @GSL_LIBS@ \
	     @BOOST_PROGOP@ @BOOST_REGEX@ @BOOST_FS@ \
	     @BOOST_SYSTEM@ @BOOST_MPI@ @BOOST_SERIALIZATION@
if NOLAUNCH
  AM_CPPFLAGS += @EO_CFLAGS@
  AM_LDFLAGS  += @EO_LIBS@
endif    

LDADD = ../opt/libopt.la 
if DOFORTRAN
   LDADD += @FCLIBS@ 
else 
if DODIRECTIAGA
   LDADD +=  ../mpi/libmpiobject.la
   AM_LDFLAGS += @FCLIBS@
endif
endif

bin_PROGRAMS = escan 

noinst_LTLIBRARIES = libescan.la
libescan_ladir = $(includedir)/LaDa/pescan_interface 
libescan_la_SOURCES = $(PESCAN_CC) 
nobase_libescan_la_HEADERS = $(PESCAN_H) $(PESCAN_IMPL_H) 
libescan_la_CPPFLAGS = $(AM_CPPFLAGS) 
if DODIRECTIAGA
  libescan_la_SOURCES += $(PESCAN_CC) $(PESCAN_F90)
  libescan_la_FCFLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  libescan_la_F77FLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  libescan_la_F90FLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
  libescan_la_FFLAGS = $(FCFLAGS) @IAGA_FCFLAGS@
endif

escan_SOURCES = main.cc
escan_CPPFLAGS = $(AM_CPPFLAGS) 
escan_LDFLAGS =  $(AM_LDFLAGS) 
escan_LDADD = libescan.la \
	      ../crystal/libcrystal.la \
	      ../atat/libatat.la \
	      ../vff/libvff.la \
              ../physics/libphysics.la \
	      ../print/libprint.la \
	      @MPI_ADDLIB@ \
	      $(LDADD) @LD_IAGA_FLAGS@



../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../crystal/libcrystal.la: 
	cd ../crystal; $(MAKE)
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
../vff/libvff.la:
	cd ../vff; $(MAKE)
../print/libprint.la:
	cd ../print; $(MAKE)
if DOFORTRAN
../mpi/libopt.la:
	cd ../opt; $(MAKE)
endif
if DODIRECTIAGA
../lapack/libmpi.la:
	cd ../mpi; $(MAKE)
endif
