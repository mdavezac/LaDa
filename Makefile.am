AUTOMAKE_OPTIONS = foreign
EXTRA_DIST = max_revision.pl
if DODARWIN 
  DARWINSUBDIR = darwin
  PESCANSUBDIR = pescan_interface
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
SUBDIRS = $(MPISUBDIR) atat physics opt lamarck print $(VFFSUBDIR) $(PESCANSUBDIR) $(DARWINSUBDIR)
revision.h: force
	./max_revision.pl
force: 
all: revision.h

