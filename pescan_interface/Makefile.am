#
# Version: $Id$
#
include ../Files
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPI_CFLAGS@ @NAG_CFLAGS@ @GSL_CFLAGS@ @BOOST_CXXFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @NAG_LIBS@ @GSL_LIBS@ @MPI_LIBS@ \
	     @BOOST_PROGOP@ @BOOST_REGEX@ @BOOST_SYSTEM@ @BOOST_MPI@
if NOLAUNCH
  AM_CPPFLAGS += @EO_CFLAGS@
  AM_LDFLAGS  += @EO_LIBS@
endif    

LDADD = 
if DOFORTRAN
   LDADD += @FCLIBS@ ../opt/libopt.la 
else 
if DODIRECTIAGA
   LDADD +=  ../mpi/libmpiobject.la
   AM_LDFLAGS += @FCLIBS@
endif
endif

if DOEMASS
 emassSOURCE = $(EMASS_CC)
 emassHEADER = $(EMASS_H)
 emassLDADD = ../lapack/liblapack.la
 binEMASS = emass
else 
 emassSOURCE = 
 emassHEADER = 
 emassLDADD = 
 bin_EMASS = 
endif

bin_PROGRAMS = escan $(binEMASS)

noinst_LTLIBRARIES = libescan.la
libescan_ladir = $(includedir)/LaDa/pescan_interface 
libescan_la_SOURCES = $(PESCAN_CC) $(emassSOURCE)
nobase_libescan_la_HEADERS = $(PESCAN_H) $(emassHEADER)
libescan_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
if DODIRECTIAGA
  libescan_la_SOURCES += $(PESCAN_CC) $(emassSOURCE) $(PESCAN_F90)
  libescan_la_FCFLAGS = $(FCFLAGS) 
endif

escan_SOURCES = main.cc
escan_CPPFLAGS = $(AM_CPPFLAGS)  -I .. 
escan_LDFLAGS =  $(AM_LDFLAGS) 
escan_LDADD = libescan.la \
	      ../crystal/libcrystal.la \
	      ../atat/libatat.la \
	      ../vff/libvff.la \
              ../physics/libphysics.la \
	      ../print/libprint.la \
	      @MPI_ADDLIB@ \
	      $(LDADD) @LD_IAGA_FLAGS@


if DOEMASS
emass_SOURCES = $(escan_SOURCES)
emass_CPPFLAGS = $(escan_CPPFLAGS) -D_EMASS
emass_LDFLAGS = $(escan_LDFLAGS)
emass_LDADD = $(escan_LDADD) $(emassLDADD)
endif

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
if DOEMASS
../lapack/liblapack.la:
	cd ../lapack; $(MAKE)
endif 
if DODIRECTIAGA
../lapack/libmpi.la:
	cd ../mpi; $(MAKE)
endif
