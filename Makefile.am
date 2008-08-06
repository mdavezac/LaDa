#
# Version: $Id$
#
AUTOMAKE_OPTIONS = foreign
EXTRA_DIST = max_revision.pl bootstrap.sh \
	     scripts/* doxy.conf docsource/* \
	     Files
include Files
if DOBANDGAPOPT 
  DARWINSUBDIR = darwin
  PESCANSUBDIR = pescan_interface
  VFFSUBDIR = vff
endif 
if DOLAYEREDOPT
  DARWINSUBDIR = darwin
  PESCANSUBDIR = pescan_interface
  VFFSUBDIR = vff
endif 
if DOEMASSOPT
  DARWINSUBDIR = darwin
  PESCANSUBDIR = pescan_interface
  VFFSUBDIR = vff
  LAPACKDIR = lapack
endif 
if DOCEOPT
  DARWINSUBDIR = darwin
endif 
if DODARWIN
  DARWINSUBDIR = darwin
endif 
if DOMPI 
  MPISUBDIR = mpi
endif 
if DOVFF
  VFFSUBDIR = vff
endif 
if DOPESCAN
  PESCANSUBDIR = pescan_interface
endif 
if DOLAPACK
  LAPACKDIR = lapack
endif
if DOSEPARABLE
  SEPARABLEDIR = separable
endif
if DOPYTHON
  PYTHONDIR = python
endif
if DOCE
  CEDIR = ce
endif
SUBDIRS = $(MPISUBDIR) atat physics opt crystal print $(CEDIR) $(VFFSUBDIR) \
	  $(PESCANSUBDIR) $(LAPACKDIR) $(DARWINSUBDIR) $(SEPARABLEDIR) \
	  $(PYTHONDIR)
revision.h: force
	@./scripts/max_revision.pl
force: ;
all: all-local
all-local: revision.h

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
endif

if DOPYTHON
python: python/liblada
python/liblada: force
	@cd python; $(MAKE) liblada.la
endif


if DODOC
doc: $(ALL_DOCS)
	@cp doxy.conf doxy.tmp
	@echo "HAVE_DOT = @FOUNDDOT@" >> doxy.tmp
	@echo "INPUT = $(ALL_DOCS)" >> doxy.tmp
	@echo "OUTPUT_DIRECTORY = @DOCDIR@ " >> doxy.tmp
	@cat doxy.tmp | doxygen -
	@rm doxy.tmp
clean-local:
	-rm -rf doc
########install-data-local: doc
########	$(INSTALL_DATA) -d $(srcdir)/doc $(docdir)
endif
