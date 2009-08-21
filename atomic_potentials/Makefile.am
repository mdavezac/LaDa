#
# Version: $Id$
#
SUBDIRS = collapse
include $(top_srcdir)/Files

AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @BOOST_CXXFLAGS@ -I $(top_srcdir) 
AM_LDFLAGS = 

noinst_LTLIBRARIES = libpotentials.la
libpotentials_la_SOURCES = $(ap_cc)
libpotentials_ladir = $(includedir)/LaDa/atomic_potentials
nobase_libpotentials_la_HEADERS = $(ap_h)
libpotentials_la_CPPFLAGS = $(AM_CPPFLAGS) 
libpotentials_la_LIBADD = collapse/libcollapse.la 

collapse/libcollapse.la: $(addprefix $(top_srcdir)/atomic_potentials/collapse/, $(ap_collapse_h) )\
                         $(addprefix $(top_srcdir)/atomic_potentials/collapse/, $(ap_collapse_cc) ) 
	@cd collapse; $(MAKE) libcollapse.la

