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
  SEPARABLEDIR = lapack
endif
SUBDIRS = $(MPISUBDIR) atat physics opt lamarck print $(VFFSUBDIR) \
	  $(PESCANSUBDIR) $(LAPACKDIR) $(DARWINSUBDIR) $(SEPARABLE)
revision.h: force
	@./scripts/max_revision.pl
force: ;
all: all-local
all-local: revision.h
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
