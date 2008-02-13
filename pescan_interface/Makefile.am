include ../Files
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPI_CFLAGS@ @NAG_CFLAGS@ @GSL_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @NAG_LIBS@ @GSL_LIBS@ @MPI_LIBS@
if NOLAUNCH
  AM_CPPFLAGS += @EO_CFLAGS@
  AM_LDFLAGS  += @EO_LIBS@
endif    

fortranLDFLAGS = 
fortranLDADD =
gslCPPFLAGS =
gslLDFLAGS = 
if DOFORTRAN
  fortranLDFLAGS =  @FCLIBS@
  fortranLDADD = ../opt/libopt.la 
endif
if DOEMASS
 fortranLDFLAGS =  @FCLIBS@
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

noinst_LTLIBRARIES = libpescan.la
libpescan_ladir = $(includedir)/pescan_interface 
libpescan_la_SOURCES = $(PESCAN_CC) $(emassSOURCE)
nobase_libpescan_la_HEADERS = $(PESCAN_H) $(emassHEADER)
libpescan_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 

escan_SOURCES = main.cc
escan_CPPFLAGS = $(AM_CPPFLAGS)  -I .. 
escan_LDFLAGS =  $(AM_LDFLAGS)
escan_LDADD = libpescan.la \
	      ../lamarck/liblamarck.la \
	      ../atat/libatat.la \
	      ../vff/libvff.la \
              ../physics/libphysics.la \
	      ../print/libprint.la \
	      @MPI_ADDLIB@ \
	      $(fortranLDADD) \
              $(fortranLDFLAGS)


if DOEMASS
emass_SOURCES = $(escan_SOURCES)
emass_CPPFLAGS = $(escan_CPPFLAGS) -D_EMASS
emass_LDFLAGS = $(escan_LDFLAGS)
emass_LDADD = $(escan_LDADD) $(emassLDADD)
endif

../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
	cd ../lamarck; $(MAKE)
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
