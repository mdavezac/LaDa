#
#  Version: $Id$
#


class Eval:

  def __init__(self, size=10):
    from random import randint
    import numpy

    self.size = size
    self.target = numpy.array([ randint(0,1] for u in xrange(self.size) ])

  def __call__(indiv):
    from math import fabs

    result = 0e0
    for r in indiv.genes - self.target:
      result += fabs(r)
    return float(result) / float(len(indiv.genes))

    
def  main():
  import lada.ga 

  darwin = object


if __name__ == "__main__":
  main()





