#
# Version: $Id$
#
#! /usr/bin/python

def allsq( _collapse, tolerance = 1e-14, maxiter=50, verbose=False ):
  """ Alternating least-square-fit function.
  """ 
  import numpy
  from random import shuffle
  from math import fabs as abs
  from lada.minimizer import cgs

  N = _collapse.nb_coordinates;
  outer_iter = 0

  old_energies = 0
  energies = 0
  while( maxiter == 0 or outer_iter < maxiter ):

    if verbose: print "allsq iteration: ", outer_iter
    outer_iter += 1
    old_energies = energies

    # Iterates over a shuffled list of directions.
    coords = range(_collapse.nb_coordinates)
    shuffle(coords)
    for i in coords:
      if verbose: print "  _ coordinate: ", i
      A, b = _collapse.lsq_data(i)
      x = _collapse.coefficients(i)
      x, res, iter = cgs(A, x, b, tolerance = tolerance)
      _collapse.update(i, x)
      if verbose: print res, iter

    
    energies = _collapse.convergence
    if verbose: print "  convergence: %18.8f %18.8f" % ( energies, energies - old_energies )
    if abs(energies - old_energies) < tolerance: break

  if verbose: print "Final convergence: %18.8f %18.8f" % ( energies, energies - old_energies )
  return outer_iter, abs(energies - old_energies)




