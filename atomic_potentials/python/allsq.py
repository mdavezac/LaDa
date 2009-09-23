#
# Version: $Id$
#
#! /usr/bin/python

def allsq( _collapse, tolerance = 1e-16, maxiter=50, verbose=False ):
  """ Alternating least-square-fit function.
  """ 
  import numpy
  from random import shuffle
  from math import fabs as abs

  N = _collapse.nb_coordinates;
  iter = 0

  old_energies = 0
  energies = 0
  while( maxiter == 0 or iter < maxiter ):

    iter += 1
    old_energies = energies

    # Iterates over a shuffled list of directions.
    coords = range(_collapse.nb_coordinates)
#   shuffle(coords)
    for i in coords:
      A, b = _collapse.lsq_data(i)
      print A
      print len(A[0]), len(A[1]), A.shape

      x = numpy.linalg.solve(A, b)
      _collapse.update(x, i)

    energies = collapse.convergence()
    if abs(energies - old_energies) < tolerance: break

  return iter, abs(energies - old_energies)




