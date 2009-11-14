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
  mlclasses = ce.MLClusterClasses("input.xml", False)
  fit = ce.PairRegulatedFit(mlclasses, alpha=3, tcoef=1)

  fit.read_directory("data")
  A, b = fit()
  x, residual, rank, s = numpy.linalg.lstsq(A, b)
  t1  = time()

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
