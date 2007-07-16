AM_CPPFLAGS = -malign-double @CPP_DEBUG_FLAGS@
AM_LDFLAGS = 
noinst_LTLIBRARIES = libmpiobject.la
libmpiobject_la_SOURCES = mpi_object.cc
libmpiobject_la_LDFLAGS = $(AM_LDFLAGS) @MPICH_LDFLAGS@
libmpiobject_la_CPPFLAGS = $(AM_CPPFLAGS) $(MPICH_CFLAGS)
