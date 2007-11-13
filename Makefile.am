AUTOMAKE_OPTIONS = foreign
EXTRA_DIST = max_revision.pl bootstrap.sh \
	     doxy.conf docsource/*
if DOVFF
  VFFSUBDIR = vff
endif 
if DOPESCAN
  PESCANSUBDIR = pescan_interface
endif 
if DOPDARWIN 
  DARWINSUBDIR = darwin
  PESCANSUBDIR = pescan_interface
  VFFSUBDIR = vff
endif 
if DOMOLECULARITY
  DARWINSUBDIR = darwin
  PESCANSUBDIR = pescan_interface
  VFFSUBDIR = vff
endif 
if DOMPI 
  MPISUBDIR = mpi
endif 
if DODARWIN
  DARWINSUBDIR = darwin
endif 
SUBDIRS = $(MPISUBDIR) atat physics opt lamarck print $(VFFSUBDIR) $(PESCANSUBDIR) $(DARWINSUBDIR)
revision.h: force
	./max_revision.pl
force: 
all: revision.h

