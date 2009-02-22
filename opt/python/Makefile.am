#
# Version: $Id$
#

if DOPYTHON
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I@PYTHON_INCLUDE_DIR@
AM_LDFLAGS = @TINYXML_LIBS@  @BOOST_LDFLAGS@ @BOOST_PYTHON@  @BOOST_FS@
if DOMPI
  AM_LDFLAGS += @BOOST_SERIALIZATION@  @BOOST_MPI@
endif

EXTRA_DIST = $(OPT_PYTHON_H) 
pkgpyexec_LTLIBRARIES = opt.la 
opt_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) 
opt_la_SOURCES = $(OPT_PYTHON_CC) 
opt_la_LDFLAGS = $(AM_LDFLAGS) -module
opt_la_LIBADD =  ../libpythonconvenience.la 
../libpythonconvenience.la:
	cd ../; $(MAKE) libpythonconvenience.la
endif
