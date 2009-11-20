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
                                B3=5, B4=5 )
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


  darwin.Individual = Individual
  # size of the individuals.
  ce.Individual.size = len(evaluation)

  darwin = Darwin()

  darwin = standard.add_population_evaluation( darwin, evaluation )
  darwin.checkpoints = [ standard.print_offspring, 
                         standard.average_fitness,
                         standard.best ]

  darwin.mating = bitstring.Mating()
  darwin.rate   = 0.2
  darwin.popsize = 40
  darwin.max_gen = 3000

  dd.run(darwin)

if __name__ == "__main__":
  main()
