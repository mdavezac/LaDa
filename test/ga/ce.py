#
#  Version: $Id$
#

def  main():
  from lada.ga import darwin as dd, bitstring, standard, ce
  from lada import crystal
  import boost.mpi
  import numpy
  import copy

  # evaluation class which keeps tracks of the number of calls.
  class Eval(ce.EvalFitPairs):
    def __init__(self, *args, **kwargs):
      ce.EvalFitPairs.__init__(self, *args, **kwargs)
      self.nbevals = 0

    def __call__(self, indiv):
      self.nbevals += 1
      return ce.EvalFitPairs.__call__(self, indiv)

  class PrintNbEvals:
    def __init__(self, evaluation): 
      self.evaluation = evaluation
    def __call__(self, darwin):
      print "Number of functional evaluations: ", self.evaluation.nbevals
      return True

  # creates  lattice
  lattice = crystal.Lattice("data/lattice.xml")


  # the evaluations (CE fitting) class.
  evaluation = Eval( lattice=lattice,
                     path = "data",
                     lmo_ratio=0.3333333, 
                     pairs = (5, 20, 5),
                     B3=9, B4=5, B5=2, B6=2)

  # the darwin class holding all ga parameters.
  class Darwin: pass
  darwin = Darwin()


  # size of the individuals.
  ce.Individual.size = len(evaluation)
  # Maximum number of bits when initializing individuals.
  ce.Individual.max_mbs_oninit = min(ce.Individual.size, 15)

  darwin.Individual = ce.Individual

  darwin = standard.add_population_evaluation( darwin, evaluation )
  darwin.checkpoints = [ standard.print_offspring, 
                         standard.average_fitness,
                         standard.best,
                         PrintNbEvals(evaluation) ]

  mating = standard.Mating(sequential=False)
  mating.add( ce.SetOver(rate=0.5, inclusive=True), rate=1 )
  mating.add( bitstring.Crossover(rate=0.5), rate=1 )
  mating.add( bitstring.Mutation(rate=2e0/float(ce.Individual.size)), rate=0.2 )

  darwin.mating = standard.Mating(sequential=True)
  darwin.mating.add( mating, rate=0.9 )
  itermax = int( float(len(evaluation)) * 1.5 )
  darwin.mating.add( bitstring.LocalSearch(evaluation, darwin, itermax=itermax), rate=0.8 )

  darwin.rate   = 0.2
  darwin.popsize = 20
  darwin.max_gen = 3000

  dd.run(darwin)

if __name__ == "__main__":
  main()
