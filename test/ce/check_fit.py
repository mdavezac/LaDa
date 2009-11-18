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

  ncls, nstr = fit.size()
  size = max( int(float(nstr)*0.3333), 1 )
  full_set = []
  predsets = []
  for i in xrange(100):
    if len(full_set) < size:
      full_set = [ j for j in range(nstr) ]
      random.shuffle(full_set)
    predsets.append( full_set[:size] )
    full_set = full_set[size:]
    print "  ", predsets[-1]

  errors = ce.leave_many_out(fit, predsets)
  prediction = 0e0
  training = 0e0
  for i in xrange(errors.shape[0]):
    for j in xrange(errors.shape[1]):
      if j in predset[i]: predictions += errors[i,j]
      else: training += errors[i,j]
  print errors.shape
  pred_x = numpy.array([j for i in s for j, s in enumerate(predsets)])
  pred_y = numpy.array([i for i in s for s in predsets] )
  prediction = errors[ pred_x, pred_y ]
  training = errors[ numpy.array([[j for i in xrange(nstr) if i not in s]\
                                  for j, s in enumerate(predsets)]),
                       numpy.array([[i for i in xrange(nstr) if i not in s] for s in predsets]) ]
  print "Training error: ", numpy.average(training*training)
  print "Prediction  error: ", numpy.average(prediction*prediction)


  errors = ce.leave_one_out(fit)
  npreds, nstr = errors.shape
  prediction = errors[ numpy.arange(errors.shape[0]), numpy.arange(errors.shape[0]) ]
  training\
    = errors[ numpy.array([[(j,i) for i in xrange(nstr) if i != j] for j in xrange(npreds)]) ]
  t2  = time()

  print x
  print "Training error: ", numpy.average(training*training)
  print "Prediction  error: ", numpy.average(prediction*prediction)
  print t2 - t1, t1 - t0

if __name__ == "__main__":
  main()
