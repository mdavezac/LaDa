#
#  Version: $Id$
#


class Eval:

  def __init__(self, size=10):
    from random import randint
    import numpy

    self.size = size
    self.target = numpy.array([randint(0,1) for u in xrange(self.size) ])

  def __call__(self, indiv):
    from math import fabs

    result = 0e0
    for r in indiv.genes - self.target:
      result += fabs(r)
    return float(result) / float(len(indiv.genes))

    
def  main():
  from lada.ga import darwin as dd, bitstring, standard
  import numpy
  import copy

  class Darwin: pass


  def stop_at_zero(self):
    for indiv in self.population:
      if indiv.fitness == 0e0: return False
    return True

  bitstring.Individual.size =30
  darwin = Darwin()

  evaluation = Eval()
  evaluation.target = numpy.array([1 for u in xrange(bitstring.Individual.size)])

  print "Target: ", evaluation.target
  darwin = dd.add_population_evaluation( darwin, evaluation )
  darwin.checkpoints = [ standard.print_offspring, 
                         standard.average_fitness,
                         standard.best,
                         stop_at_zero ]
  darwin.mating = bitstring.Mating()
  darwin.rate   = 0.2
  darwin.popsize = 100
  darwin.max_gen = 300

  dd.run(darwin)

if __name__ == "__main__":
  main()
