#
#  Version: $Id$
#

def main():

  from boost import mpi
  from lada import crystal, ce
  import numpy
  from math import fabs
  from time import time


  t0  = time()
  
  lattice = crystal.Lattice("input.xml")
  mlclasses = ce.create_clusters( lattice, nth_shell=0, order=0, site=0)
# mlclasses.extend( ce.create_clusters( lattice, nth_shell=0, order=1, site=0) )
  mlclasses.extend( ce.create_clusters( lattice, nth_shell=20, order=2, site=0) )
  mlclasses.extend( ce.create_clusters( lattice, nth_shell=9, order=3, site=0) )
  mlclasses.extend( ce.create_clusters( lattice, nth_shell=4, order=4, site=0) )
  mlclasses.extend( ce.create_clusters( lattice, nth_shell=2, order=5, site=0) )
  mlclasses.extend( ce.create_clusters( lattice, nth_shell=2, order=6, site=0) )
  print mlclasses, len(mlclasses)

  t1 = time()
  print t1 - t0

if __name__ == "__main__":
  main()
