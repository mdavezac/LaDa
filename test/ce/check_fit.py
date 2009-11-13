#
#  Version: $Id$
#


def main():

  from boost import mpi
  from lada import crystal, ce
  import numpy.linalg


  lattice = crystal.Lattice("input.xml")
  mlclasses = ce.MLClusterClasses("input.xml", False)
  fit = ce.PairRegulatedFit(mlclasses)

  fit.read_directory("data")

  print fit._pis.shape
  print len(fit._str_onoff)
  print len(fit._cls_onoff)


  A,b = fit()
  x, residues, rank, s = numpy.linalg.lstsq(A, b)
  print x
  print residues

if __name__ == "__main__":
  main()
