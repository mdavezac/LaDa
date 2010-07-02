#
# Version: $Id$
#
ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign
include $(top_srcdir)/Files
EXTRA_DIST = $(top_srcdir)/bootstrap.sh \
	     $(EXTRA_FILES) \
             $(top_srcdir)/doxy.conf \
	     $(top_srcdir)/Files
if DOMPI 
  MPISUBDIR = mpi
endif 
if DOVFF
  VFFSUBDIR = vff
endif 
if DOPESCAN
  PESCANSUBDIR = escan
endif 
if DOSEPARABLE
  SEPARABLEDIR = separable
  CEDIR = ce
endif
if DOPYTHON
  PYTHONDIR = python ga jobs
if DOVASP
  PYTHONDIR += vasp
endif 
endif
if DOCE
  CEDIR = ce
endif
if DOMODELS
  MODELSDIR = models
endif
if DOPOTENTIALS
  POTENTIALSDIR = atomic_potentials
endif
if DOENUM
  ENUMDIR = enum
endif
SUBDIRS = $(PYTHONDIR) $(MPISUBDIR) math physics opt minimizer crystal print $(CEDIR) $(VFFSUBDIR) \
	  $(PESCANSUBDIR) $(SEPARABLEDIR) $(MODELSDIR) $(POTENTIALSDIR) $(ENUMDIR) 
# revision.h: force
# 	@cd $(top_srcdir); scripts/max_revision.pl
# force: ;
# all: all-local
# all-local: revision.h

if DOCE
ce/cubic: force
	@cd ce; $(MAKE) cubic
ce/tetra: force
	@cd ce; $(MAKE) tetra
ce/regular: force
	@cd ce; $(MAKE) regular
endif

if DOSEPARABLE
separable/cubic: force
	@cd separable; $(MAKE) cubic
separable/tetra: force
	@cd separable; $(MAKE) tetra
separable/fx/cubic: force
	@cd separable/fixedlattice; $(MAKE) cubic
separable/fx/tetra: force
	@cd separable/fixedlattice; $(MAKE) tetra
separable/fx/mixed: force
	@cd separable/fixedlattice; $(MAKE) mixed
separable/fx/mixedt: force
	@cd separable/fixedlattice; $(MAKE) mixedtetra
endif
if DOPYTHON 
python: force
	@cd python; $(MAKE)
endif

if DOPYTHON
python: python/liblada
python/liblada: force
	@cd python; $(MAKE) liblada.la
endif


if DODOC
doc: $(ALL_DOCS)
	@cp $(top_srcdir)/doxy.conf doxy.tmp
	@echo "HAVE_DOT = @FOUNDDOT@" >> doxy.tmp
	@echo "INPUT = $(ALL_DOCS)" >> doxy.tmp
	@cat doxy.tmp | doxygen -
	@rm doxy.tmp
clean-local:
	-rm -rf doc
########install-data-local: doc
########	$(INSTALL_DATA) -d $(srcdir)/doc $(docdir)
endif

