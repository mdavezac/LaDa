#
# Version: $Id$
#
#! /usr/bin/python

def allsq( _collapse, tolerance = 1e-14, itermax=50, verbose=False, cgsitermax=10, cgstolerance=1e-14 ):
  """ Alternating least-square-fit function.
  """ 
  import numpy
  from random import shuffle
  from math import fabs as abs, sqrt
  from scipy.sparse.linalg import cgs

  N = _collapse.nb_coordinates;
  outer_iter = 0

  old = 0
  old_energies = 0
  energies = _collapse.convergence
  while( itermax == 0 or outer_iter < itermax ):

    if verbose: print "allsq iteration: ", outer_iter
    outer_iter += 1
    old_energies = energies

    # Iterates over a shuffled list of directions.
    coords = range(_collapse.nb_coordinates)
    shuffle(coords)
    for i in coords:
      if verbose: print "  _ coordinate: ", i,
      A, b = _collapse.lsq_data(i)
      Ac, b = _collapse.lsq_data(i)
      x = _collapse.coefficients(i)
      if verbose:  old = _collapse.convergence
      x, info = cgs(A, b, x, tol=cgstolerance, maxiter=cgsitermax)
      _collapse.update(i, x)
      if verbose:
        n =  _collapse.convergence
        print info, n
    #   assert old > n or abs(old - n) < 1e-12,  "Fit is degrading: %e < %e.\n" % (old, n)

    

    energies = _collapse.convergence
    if verbose: print "  convergence: %18.8f %18.8f" % ( energies, energies - old_energies )
    else: print "  convergence: %18.8f %18.8f" % ( energies, energies - old_energies )
   #assert old_energies >  energies or abs(old_energies -  energies) < 1e-12, \
   #       "Fit is degrading: %e < %e.\n" % (old_energies, _energies)
    if abs(energies - old_energies) < tolerance: break

  if verbose: print "Final convergence: %18.8f %18.8f" % ( energies, energies - old_energies )
  return outer_iter, energies




