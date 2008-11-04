#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPI_CFLAGS@ @EO_CFLAGS@ -I $(top_srcdir) @BOOST_CXXFLAGS@ 
AM_LDFLAGS = @TINYXML_LIBS@ @EO_LIBS@ @BOOST_LDFLAGS@ @BOOST_MPI@ \
	     @BOOST_SERIALIZATION@ @BOOST_PROGOP@ @BOOST_FS@
if DOBANDGAPOPT
  BANDGAPOPT=bandgap_opt
endif
if DOCEOPT
  CEOPT = ce_opt
endif
if DOLAYEREDOPT
  LAYEREDOPT = layered_opt
endif
if DOALLOYLAYERSOPT
  ALLOYLAYERSOPT = alloylayers_opt
endif
LDADD =
if DOFORTRAN
  LDADD +=  @FCLIBS@
endif
if DODIRECTIAGA
   AM_LDFLAGS += @FCLIBS@
   LDADD += @FCLIBS@
endif
bin_PROGRAMS = $(CEOPT) $(BANDGAPOPT) $(LAYEREDOPT) $(ALLOYLAYERSOPT)
ce_opt_SOURCES = main.cc $(DARWIN_CE_CC)
ce_opt_HEADERS = $(DARWIN_CE_H)
ce_optdir = $(includedir)/LaDa/darwin
ce_opt_CPPFLAGS =  -D _CE $(AM_CPPFLAGS) @GSL_CFLAGS@ -D _CUBIC_CE_
ce_opt_LDFLAGS =  $(AM_LDFLAGS)
ce_opt_LDADD = libdarwin.la \
	       ../ce/libce.la \
	       ../crystal/libcrystal.la \
	       ../physics/libphysics.la \
	       ../atat/libatat.la \
	       ../print/libprint.la \
	       ../mpi/libmpiobject.la \
	       ../opt/libopt.la \
	       @MPI_ADDLIB@

bandgap_opt_SOURCES = main.cc $(DARWIN_BG_CC)
bandgap_opt_HEADERS = $(DARWIN_BG_H) 
bandgap_optdir = $(includedir)/LaDa/darwin
bandgap_opt_CPPFLAGS =  -D _PESCAN $(AM_CPPFLAGS) @NAG_CFLAGS@ @GSL_CFLAGS@ 
bandgap_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIBS@ @NAG_LIBS@ 
bandgap_opt_LDADD = $(ce_opt_LDADD) \
                    ../vff/libvff.la \
	            ../pescan_interface/libescan.la \
		    ../physics/libphysics.la \
		    $(fortranLDADD) \
	            @MPI_ADDLIB@ \
	            $(LDADD) @LD_IAGA_FLAGS@

layered_opt_SOURCES = main.cc $(DARWIN_L_CC)
layered_opt_HEADERS = $(DARWIN_L_H) 
layered_optdir = $(includedir)/darwin
layered_opt_CPPFLAGS =  -D _MOLECULARITY $(AM_CPPFLAGS) @NAG_CFLAGS@ @GSL_CFLAGS@ 
layered_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIBS@ @NAG_LIBS@ 
layered_opt_LDADD = $(bandgap_opt_LDADD)

alloylayers_opt_SOURCES = main.cc $(DARWIN_ALLOYLAYERS_CC)
alloylayers_opt_HEADERS = $(DARWIN_ALLOYLAYERS_H) 
alloylayers_optdir = $(includedir)/darwin
alloylayers_opt_CPPFLAGS =  -D _ALLOY_LAYERS_ $(AM_CPPFLAGS) @NAG_CFLAGS@ @GSL_CFLAGS@ 
alloylayers_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIBS@ @NAG_LIBS@ @BOOST_SIGNAL@
alloylayers_opt_LDADD = $(bandgap_opt_LDADD)

noinst_LTLIBRARIES = libdarwin.la
libdarwin_ladir = $(includedir)/LaDa/darwin
libdarwin_la_SOURCES = $(DARWIN_LIB_CC)
nobase_libdarwin_la_HEADERS = $(DARWIN_LIB_H) 
libdarwin_la_CPPFLAGS = $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../ce/libce.la: 
	cd ../ce; $(MAKE)
../crystal/libcrystal.la: 
	cd ../crystal; $(MAKE)
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
../pescan_interface/libescan.la:
	cd ../pescan; $(MAKE)
../print/libprint.la:
	cd ../print; $(MAKE)
../opt/libopt.la:
	cd ../opt; $(MAKE)
