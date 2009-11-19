#
#  Version: $Id$
#


def main():

  from boost import mpi
  from lada import crystal, ce
  import numpy
  from math import fabs
  from time import time
  import random


  t0  = time()
  lattice = crystal.Lattice("input.xml")
  mlclasses = ce.MLClusterClasses("input.xml", False)
  fit = ce.PairRegulatedFit(mlclasses, alpha=3, tcoef=1)

  fit.read_directory("data")
  A, b = fit()
  x, residual, rank, s = numpy.linalg.lstsq(A, b)
  t1  = time()

  t2  = time()
  loo_errors = ce.leave_one_out(fit)
  npreds, nstr = loo_errors.shape
  prediction = 0e0
  training = 0e0
  for i in xrange(loo_errors.shape[0]):
    for j in xrange(loo_errors.shape[1]):
      if j == i: prediction += loo_errors[i,j]*loo_errors[i,j]
      else: training += loo_errors[i,j]*loo_errors[i,j]
  prediction /= float(nstr)
  training /= float(nstr*nstr - nstr)
  print "Training error: ", training
  print "Prediction  error: ", prediction

  ncls, nstr = fit.size()
  size = max( int(float(nstr)*0.33333333333), 1 )
  full_set = []
  predsets = []
  for i in xrange(10):
    if len(full_set) < size:
      full_set = [ j for j in range(nstr) ]
      random.shuffle(full_set)
    predsets.append( full_set[:size] )
    full_set = full_set[size:]

  errors = ce.leave_many_out(fit, predsets)
  prediction = 0e0
  training = 0e0
  n = 0
  for i in xrange(errors.shape[0]):
    for j in xrange(errors.shape[1]):
      if j in predsets[i]:
        prediction += errors[i,j]*errors[i,j]
        n += 1
      else: training += errors[i,j]*errors[i,j]
  prediction /= float(n)
  training /= float(errors.shape[0]*errors.shape[1] - n)
  print errors.shape, n
  print "Training error: ", training
  print "Prediction  error: ", prediction



  print t2 - t1, t1 - t0

if __name__ == "__main__":
  main()
