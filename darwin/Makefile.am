#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @MPI_CFLAGS@ @EO_CFLAGS@ -I $(top_srcdir) @BOOST_CXXFLAGS@ 
AM_LDFLAGS = @TINYXML_LIBS@ @EO_LIBS@ @BOOST_LDFLAGS@ @BOOST_MPI@ \
	     @BOOST_SERIALIZATION@ @BOOST_PROGOP@ @BOOST_FS@ @BOOST_SIGNAL@
EXTRA_DIST = $(top_srcdir)/darwin/main.extras.h
if DOCEOPT
  CEOPT = ce_opt
endif
if DOALLOYLAYERSOPT
  ALLOYLAYERSOPT = alloylayers_opt
endif
LDADD = @FCLIBS@
if DODIRECTIAGA
   AM_LDFLAGS += @FCLIBS@
endif
bin_PROGRAMS = $(CEOPT) $(ALLOYLAYERSOPT)
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
		../minimizer/libminimizer.la \
	       @MPI_ADDLIB@

alloylayers_opt_SOURCES = main.cc 
alloylayers_opt_HEADERS = alloylayers/main.extras.h
alloylayers_optdir = $(includedir)/LaDa/darwin
alloylayers_opt_CPPFLAGS =  -D _ALLOY_LAYERS_ $(AM_CPPFLAGS) @GSL_CFLAGS@ 
alloylayers_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIBS@ @BOOST_SIGNAL@
alloylayers_opt_LDADD = libdarwin.la \
			liballoylayers.la \
			../vff/libvff.la \
			../pescan_interface/libescan.la \
	                ../crystal/libcrystal.la \
	                ../physics/libphysics.la \
	                ../atat/libatat.la \
	                ../print/libprint.la \
	                ../mpi/libmpiobject.la \
	                ../opt/libopt.la \
		        ../minimizer/libminimizer.la \
	                @MPI_ADDLIB@ \
	                $(LDADD) @LD_IAGA_FLAGS@

noinst_LTLIBRARIES = libdarwin.la libdarwinoperator.la liballoylayers.la \
                     libdarwingraph.la libdarwinckp.la libdarwintaboos.la

libdarwinoperator_ladir = $(includedir)/LaDa/darwin/operators
libdarwinoperator_la_SOURCES = $(DARWIN_OP_CC)
libdarwinoperator_la_HEADERS = $(DARWIN_OP_H)
libdarwinoperator_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwingraph_ladir = $(includedir)/LaDa/darwin/graph
libdarwingraph_la_SOURCES = $(DARWIN_GRAPH_CC)
libdarwingraph_la_HEADERS = $(DARWIN_GRAPH_H)
libdarwingraph_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwinckp_ladir = $(includedir)/LaDa/darwin/checkpoints
libdarwinckp_la_SOURCES = $(DARWIN_CKP_CC)
libdarwinckp_la_HEADERS = $(DARWIN_CKP_H)
libdarwinckp_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwintaboos_ladir = $(includedir)/LaDa/darwin/taboos
libdarwintaboos_la_SOURCES =
libdarwintaboos_la_HEADERS = $(DARWIN_TABOOS_H)
libdarwintaboos_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwinbitstring_ladir = $(includedir)/LaDa/darwin/bitstring
libdarwinbitstring_la_HEADERS = $(DARWIN_BITSTRING_H)
libdarwinbitstring_la_CPPFLAGS = $(AM_CPPFLAGS)

liballoylayers_ladir = $(includedir)/LaDa/darwin/alloylayers
liballoylayers_la_SOURCES = $(LIB_ALLOYLAYERS_CC)
liballoylayers_la_HEADERS = $(LIB_ALLOYLAYERS_H)
liballoylayers_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwin_ladir = $(includedir)/LaDa/darwin
libdarwin_la_SOURCES = $(DARWIN_LIB_CC)
nobase_libdarwin_la_HEADERS = $(DARWIN_LIB_H) 
libdarwin_la_CPPFLAGS = $(AM_CPPFLAGS) 
libdarwin_la_LIBADD = libdarwinoperator.la libdarwingraph.la libdarwinckp.la

../atat/libatat.la: 
	cd ../atat; $(MAKE) libatat.la
../ce/libce.la: 
	cd ../ce; $(MAKE) libce.la
../crystal/libcrystal.la: 
	cd ../crystal; $(MAKE) libcrystal.la
../physics/libphysics.la:
	cd ../physics; $(MAKE) libphysics.la
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE) libmpi.la
../pescan_interface/libescan.la:
	cd ../pescan; $(MAKE) libpescan.la
../print/libprint.la:
	cd ../print; $(MAKE) libprint.la
../opt/libopt.la:
	cd ../opt; $(MAKE) libopt.la
../minimizer/libminimizer.la:
	cd ../minimizer; $(MAKE) libminimizer.la
