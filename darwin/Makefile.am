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
ce_opt_SOURCES = main.cc $(darwin_ce_cc)
ce_opt_HEADERS = $(darwin_ce_h)
ce_optdir = $(includedir)/LaDa/darwin
ce_opt_CPPFLAGS =  -D _CE $(AM_CPPFLAGS) @GSL_CFLAGS@ -D _CUBIC_CE_
ce_opt_LDFLAGS =  $(AM_LDFLAGS)
ce_opt_LDADD = libdarwin.la \
	       $(top_builddir)/ce/libce.la \
	       $(top_builddir)/crystal/libcrystal.la \
	       $(top_builddir)/physics/libphysics.la \
	       $(top_builddir)/atat/libatat.la \
	       $(top_builddir)/print/libprint.la \
	       $(top_builddir)/mpi/libmpiobject.la \
	       $(top_builddir)/opt/libopt.la \
	       $(top_builddir)/minimizer/libminimizer.la \
	       @MPI_ADDLIB@

alloylayers_opt_SOURCES = main.cc 
alloylayers_opt_HEADERS = alloylayers/main.extras.h
alloylayers_optdir = $(includedir)/LaDa/darwin
alloylayers_opt_CPPFLAGS =  -D _ALLOY_LAYERS_ $(AM_CPPFLAGS) @GSL_CFLAGS@ 
alloylayers_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIBS@ @BOOST_SIGNAL@
alloylayers_opt_LDADD = libdarwin.la \
			liballoylayers.la \
			$(top_builddir)/vff/libvff.la \
			$(top_builddir)/pescan_interface/libescan.la \
	                $(top_builddir)/crystal/libcrystal.la \
	                $(top_builddir)/physics/libphysics.la \
	                $(top_builddir)/atat/libatat.la \
	                $(top_builddir)/print/libprint.la \
	                $(top_builddir)/mpi/libmpiobject.la \
	                $(top_builddir)/opt/libopt.la \
		        $(top_builddir)/minimizer/libminimizer.la \
	                @MPI_ADDLIB@ \
	                $(LDADD) @LD_IAGA_FLAGS@

noinst_LTLIBRARIES = libdarwin.la libdarwinoperator.la liballoylayers.la \
                     libdarwingraph.la libdarwinckp.la libdarwintaboos.la

libdarwinoperator_ladir = $(includedir)/LaDa/darwin/operators
libdarwinoperator_la_SOURCES = $(darwin_op_cc)
libdarwinoperator_la_HEADERS = $(darwin_op_h) $(darwin_op_impl_h)
libdarwinoperator_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwingraph_ladir = $(includedir)/LaDa/darwin/graph
libdarwingraph_la_SOURCES = $(darwin_graph_cc)
libdarwingraph_la_HEADERS = $(darwin_graph_h)
libdarwingraph_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwinckp_ladir = $(includedir)/LaDa/darwin/checkpoints
libdarwinckp_la_SOURCES = $(darwin_ckp_cc)
libdarwinckp_la_HEADERS = $(darwin_ckp_h)
libdarwinckp_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwintaboos_ladir = $(includedir)/LaDa/darwin/taboos
libdarwintaboos_la_SOURCES =
libdarwintaboos_la_HEADERS = $(darwin_taboos_h)
libdarwintaboos_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwinbitstring_ladir = $(includedir)/LaDa/darwin/bitstring
libdarwinbitstring_la_HEADERS = $(darwin_bitstring_h)
libdarwinbitstring_la_CPPFLAGS = $(AM_CPPFLAGS)

liballoylayers_ladir = $(includedir)/LaDa/darwin/alloylayers
liballoylayers_la_SOURCES = $(lib_alloylayers_cc)
liballoylayers_la_HEADERS = $(LIB_ALLOYLAYERS_H)
liballoylayers_la_CPPFLAGS = $(AM_CPPFLAGS)

libdarwin_ladir = $(includedir)/LaDa/darwin
libdarwin_la_SOURCES = $(darwin_lib_cc)
nobase_libdarwin_la_HEADERS = $(darwin_lib_h) 
libdarwin_la_CPPFLAGS = $(AM_CPPFLAGS) 
libdarwin_la_LIBADD = libdarwinoperator.la libdarwingraph.la libdarwinckp.la
