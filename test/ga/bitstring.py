#
#  Version: $Id$
#


class Eval:

  def __init__(self, size=10):
    from random import randint
    import numpy

    self.size = size
    self.target = numpy.array([randint(0,1) for u in xrange(self.size) ])
    self.nbevals = 0

  def __call__(self, indiv):
    from math import fabs

    self.nbevals += 1
    result = 0e0
    for r in indiv.genes - self.target:
      result += fabs(r)
    return float(result) / float(len(indiv.genes))

class PrintNbEvals:
  def __init__(self, evaluation): 
    self.evaluation = evaluation
  def __call__(self, darwin):
    print "Number of functional evaluations: ", self.evaluation.nbevals
    return True


    
def  main():
  from lada.ga import darwin as dd, bitstring, standard
  import numpy
  import copy

  class Darwin: pass


  def stop_at_zero(self):
    for indiv in self.population:
      if indiv.fitness == 0e0: return False
    return True

  bitstring.Individual.size =50
  darwin = Darwin()

  evaluation = Eval()
  evaluation.target = numpy.array([1 for u in xrange(bitstring.Individual.size)])
  print "Target: ", evaluation.target

  darwin = standard.add_population_evaluation( darwin, evaluation )
  darwin.checkpoints = [ standard.print_offspring, 
                         standard.average_fitness,
                         standard.best,
                         PrintNbEvals(evaluation),
                         stop_at_zero ]

  mating = standard.Mating(sequential=False)
  mating.add( bitstring.Crossover(rate=0.5), rate=0.8 )
  mating.add( bitstring.Mutation(rate=3e0/float(bitstring.Individual.size)), rate=0.2 )

  darwin.mating = standard.Mating(sequential=True)
  darwin.mating.add( mating, rate=0.8 )
  darwin.mating.add( bitstring.LocalSearch(evaluation, darwin, itermax=75), rate=0.8 )

  darwin.rate   = 0.2
  darwin.popsize = 100
  darwin.max_gen = 300

  dd.run(darwin)

if __name__ == "__main__":
  main()
