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
  evaluation = ce.Eval( lattice=lattice,
                        path = "data",
                        lmo_ratio=0.3333333, 
                        alpha=2.0, tcoef=10,
                        2B=20, 3B9, 4B=5, 5B=2, 6B=2)

  # the darwin class holding all ga parameters.
  class Darwin: pass

  # individual type.
  darwin.Individual = ce.Individual
  # size of the individuals.
  darwin.Individual.size = len(evaluation)

  darwin = Darwin()

  darwin = standard.add_population_evaluation( darwin, evaluation )
  darwin.checkpoints = [ standard.print_offspring, 
                         standard.average_fitness,
                         standard.best ]

  darwin.mating = bitstring.Mating()
  darwin.rate   = 0.2
  darwin.popsize = 100
  darwin.max_gen = 3000

  dd.run(darwin)

if __name__ == "__main__":
  main()
