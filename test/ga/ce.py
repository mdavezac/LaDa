#
#  Version: $Id$
#

def  main():
  from lada.ga import darwin as dd, bitstring, standard, ce
  from lada import crystal
  import numpy
  import copy

  # creates  lattice
  lattice = crystal.Lattice("data/lattice.xml")

  # the evaluations (CE fitting) class.
  evaluation = ce.EvalFitPairs( lattice=lattice,
                                path = "data",
                                lmo_ratio=0.3333333, 
                                pairs = (5, 20, 5),
                                B3=9, B4=5, B5=2, B6=2)
                                # B2=10) # , B3=9, B4=5, B5=2, B6=2)

  # the darwin class holding all ga parameters.
  class Darwin: pass
  darwin = Darwin

  # individual type.
  class Individual(ce.Individual):

    def __init__(self):
      from random import shuffle, randint

      ce.Individual.__init__(self)

      list_ = [i for i in xrange(len(self.genes))]
      shuffle(list_)
      self.genes[:] = False
      self.genes[list_[:randint(1, min(20, ce.Individual.size))]] = True

  class Mating(bitstring.Mating): 
    """ Adds a different crossover to the bitstring crossover """
    def __init__(self, setover = 0.8, **kwargs):

      bitstring.Mating.__init__(self, **kwargs)

      class SetOver:
        def __init__(self, rate=0.5, inclusive=True): self.rate = rate; self.inclusive = inclusive
        def __call__(self, a, b): 
          from copy import deepcopy
          from random import random

          result = deepcopy(a)
          set_a = set()
          for i, n in enumerate(a.genes):
            if n: set_a.add(i)
          set_b = set()
          for i, n in enumerate(b.genes):
            if n: set_b.add(i)

          if self.inclusive: 
            result.genes[:] = False
            for i in set_a & set_b: 
              result.genes[i] = True
            for i in set_a ^ set_b: 
              if random() < self.rate:
                result.genes[i] = True
          else: 
            for i in set_a | set_b: 
              if random() < self.rate:
                result.genes[i] = True
          return result;

      norm = 1 + setover
      self.mutrate = (1-self.rate)/norm
      self.crossrate = self.rate/norm
      self.setrate = setover/norm

      self.setover = SetOver()

    def __call__(self, darwin):

      from random import random
      
      a = darwin.selection(darwin)
      
      indiv = None
      r = random() 
      if r < self.mutrate: indiv = self.mutation(darwin.population[a])
      elif r < self.crossrate:
        b = a
        while( b == a ): b = darwin.selection(darwin)
        indiv = self.crossover( darwin.population[a], darwin.population[b])
      else: 
        b = a
        while( b == a ): b = darwin.selection(darwin)
        indiv = self.setover( darwin.population[a], darwin.population[b])
      
      if hasattr(indiv, "fitness"): delattr(indiv, "fitness")
      return indiv



  darwin.Individual = Individual
  # size of the individuals.
  ce.Individual.size = len(evaluation)

  darwin = Darwin()

  darwin = standard.add_population_evaluation( darwin, evaluation )
  darwin.checkpoints = [ standard.print_offspring, 
                         standard.average_fitness,
                         standard.best ]

  darwin.mating = Mating(rate=0.5, setover=3)
  darwin.rate   = 0.2
  darwin.popsize = 20
  darwin.max_gen = 3000

  dd.run(darwin)

if __name__ == "__main__":
  main()
