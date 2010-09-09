""" Genetic Algorithm Package
    
    The package consists of three main modules, darwin, standard, and
    bitstring, as well as application modules.  It is articulated around a user
    defined class (Darwin throughout) and the algorithm subroutine
    darwin.run(...). In general, the user will only need add to the Darwin an
    evaluation function, an Individual class type, and mating operators. Other
    components, such as selection rules for parents, taboo operators,
    population size, checkpoints and so forth are added automatically by
    darwin.run() using standard operations defined in standard using
    standard.fill_attributes(...).  
    
    The following describes the standard operations in a self instance of
    Darwin. Note that at present, none of the methods added to Darwin are bound
    to the instance. E.g. they take self explicitely:
      - Selection is performed using the self.selection method. It takes one
        argument, self itself.
      - self.ating operations take the Darwin instance as an argument. It should
        return an instance of self.Individual for which no fitness attribute
        exist. Examples are given in the bitstring module.
      - self.checkpoints is a list of unbound method taking self as the
        only argument and returning False if the GA should stop. They are executed
        at the end of each generation. They are used for both output (always
        returning True) and convergence. A method checking for maximum number
        of generations (self.max_gen) is ALWAYS added to these checkpoints.
      - self.popsize is the population size.
      - self.rate * self.popsize is the number of offspring at each
        generation. At this point, all offspring are indiscreminately entered in
        the next population.
      - self.cmp_indiv compares the fitness of two individuals. By default,
        this GA minimizes.
      - self.Individual is a class (not an instance) defining individuals. It
        should contain an __init__ function which randomly initializes the
        individual. By default, and instance of self.Individual should have
        an fitness attribute if and only if it has already been evaluated.
      - self.evaluation is a function taking self as its only argument. It
        should evaluate the offspring and populations as necessary. When using
        standard.add_population_evaluation(self, evaluation), each individual
        of population and offspring is evaluated only if that individual does
        not possess a fitness attribute. 

    An example (function test()) can be found directly in the package's __init__.py
"""

__all__ = [@which_gapackages@]

def test():
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

  darwin = standard.add_population_evaluation( darwin, evaluation )
  darwin.checkpoints = [ standard.print_offspring, 
                         standard.average_fitness,
                         standard.best,
                         stop_at_zero ]
  darwin.mating = bitstring.Mating()
  darwin.rate   = 0.2
  darwin.popsize = 100
  darwin.max_gen = 300

  dd.run(darwin)
       

if __name__ == "__main__": test()
