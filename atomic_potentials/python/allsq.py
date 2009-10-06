#
# Version: $Id$
#
#! /usr/bin/python

class Allsq(object):
  """ Alternating least-square-fit function.
  """ 

  def __init__( self, tolerance = 1e-14, itermax=50, verbose=False, \
                cgsitermax=10, cgstolerance=1e-14, bestof = 1 ):

    self.tolerance = tolerance 
    self.itermax = itermax 
    self.verbose = verbose 
    self.cgsitermax = cgsitermax
    self.cgstolerance = cgstolerance
    self.bestof = bestof

  def __call__(self, _collapse, **kwargs):
    from math import fabs
    import numpy
    import random

    if len(kwargs) != 0:
      other = Allsq(self.tolerance, self.itermax, self.verbose,\
                    self.cgsitermax, self.cgstolerance, self.bestof)
      for key in kwargs:
        if not hasattr(other, key ):
          raise AttributeError, "Unknown keyword: %s=%s.\n" % (key, kwargs[key])
        other.__setattr__(key, kwargs[key])
      return other(_collapse)

    if self.bestof < 2: return self.call(_collapse)

    best_coefs, best_s, best_e, best_iter, which = [], [], None, None, 0
    for i in range(self.bestof):
      # randomize sumofseps
      for x in _collapse.scales: x=1e0
      for j in  range(_collapse.nb_coordinates):
        x = _collapse.coefficients(j)
        for u in x: u = random.uniform(0.9, 1.1)
        _collapse.update(j, x) 

      iter, e = self.call(_collapse)
      if best_e == None or best_e > e:
        which = i
        best_e = e
        best_iter = iter
        best_coefs = []
        best_s = [x for x in _collapse.scales]
        for j in range(_collapse.nb_coordinates):
          best_coefs.append( _collapse.coefficients(j).copy() )

    if which != self.bestof:
      for i,x in enumerate(_collapse.scales): x=best_s[i]
      for i, coefs in enumerate(best_coefs):
        _collapse.update(i, coefs)

    assert fabs( _collapse.convergence - best_e ) < 1e-8, "%f %f" % (_collapse.convergence, best_e)
    return best_iter, best_e

  def call(self, _collapse):
    import numpy
    from random import shuffle
    from math import fabs as abs, sqrt
    from scipy.sparse.linalg import cgs
    from scipy.linalg import lstsq
    
    N = _collapse.nb_coordinates;
    outer_iter = 0
    
    old = 0
    old_energies = 0
    energies = _collapse.convergence
    while( self.itermax == 0 or outer_iter < self.itermax ):
    
      if self.verbose: print "allsq iteration: ", outer_iter
      outer_iter += 1
      old_energies = energies
    
      # Iterates over a shuffled list of directions.
      coords = range(_collapse.nb_coordinates)
      shuffle(coords)
      for i in coords:
        if self.verbose: print "  _ coordinate: ", i,
        A, b = _collapse.lstsq_data(i)
        x = _collapse.coefficients(i)
        # x, info = cgs(A, b, x, tol=self.cgstolerance, maxiter=self.cgsitermax)
        x, res, rank, svd = lstsq(A, b, cond=self.cgstolerance, overwrite_a=1, overwrite_b=1)
        _collapse.update(i, x)
        if self.verbose:
          n =  _collapse.convergence
          print res, n
      #   assert old > n or abs(old - n) < 1e-12,  "Fit is degrading: %e < %e.\n" % (old, n)
    
      
    
      energies = _collapse.convergence
      if self.verbose:
        print "  convergence: %18.8f %18.8f" % ( energies, energies - old_energies )
     #else: print "  convergence: %18.8f %18.8f" % ( energies, energies - old_energies )
     #assert old_energies >  energies or abs(old_energies -  energies) < 1e-12, \
     #       "Fit is degrading: %e < %e.\n" % (old_energies, _energies)
      if abs(energies - old_energies) < self.tolerance: break
    
    if self.verbose:
      print "Final convergence: %18.8f %18.8f" % ( energies, energies - old_energies )
    return outer_iter, energies


