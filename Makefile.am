AUTOMAKE_OPTIONS = foreign
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
SUBDIRS = $(MPISUBDIR) atat physics opt lamarck  $(VFFSUBDIR) $(PESCANSUBDIR) $(DARWINSUBDIR)
