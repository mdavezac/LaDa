AUTOMAKE_OPTIONS = foreign
EXTRA_DIST = max_revision.pl bootstrap.sh \
	     scripts/* doxy.conf docsource/*
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
SUBDIRS = $(MPISUBDIR) atat physics opt lamarck print $(VFFSUBDIR) \
	  $(PESCANSUBDIR) $(LAPACKDIR) $(DARWINSUBDIR) 
revision.h: force
	@scripts/max_revision.pl
force: 
all: revision.h

