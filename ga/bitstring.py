#
#  Version: $Id$
#

class Individual:
  """ An individual for bitstrings.
      The size of the individual is given statically by Individual.size
  """

  size = 10

  def __init__(self):
    """ Initializes a bitstring individual randomly.
    """
    from random import randint
    import numpy

    self.genes = numpy.array([ randint(0,1) for i in xrange(Individual.size) ])
  
  def __eq__(self, a): 
    from math import fabs
    if a == None: return False
    for i, u in enumerate(a.genes):
      if fabs( self.genes[i] - u ) > 1e-18: return False
    return True

  def __str__(self): return "%s" % (self.genes)


class Crossover:
  """ A crossover operation.
  """

  def __init__(self, rate = 0.5):
    self.rate = rate

  def __call__(self, a, b):
    from copy import deepcopy
    
    at = self.rate * len(a.genes) 
    a.genes[at:] = b.genes[len(b.genes)-at:]

    if hasattr(a, "fitness"): delattr(a, "fitness")
    return a

class Mutation:
  """ A crossover operation.
  """

  def __init__(self, rate = 0.10):
    self.rate = rate

  def __call__(self, a):
    from copy import deepcopy
    from random import uniform
    
    for i in xrange(len(a.genes)):
      if uniform(0, 1)  < self.rate:
        if a.genes[i] == 1: a.genes[i] = 0
        else: a.genes[i] = 1
    if hasattr(a, "fitness"): delattr(a, "fitness")
    return a


class LocalSearch(object):
  """ Performs a local search over a bitstring by flipping random bits. """

  def __init__(self, evaluation, darwin, itermax=3, decrease=False):
    """ Initializes a LocalSearch instance.
          _ evaluation is a functor or function taking an individual as its argument. 
          _ darwin is a class containing a taboo, a selection, and an cmp_indiv
            procedure, as well as a population.
          _ itermax is the maximum number of evaluation. 
    """

    self.evaluation = evaluation
    self.darwin = darwin
    self.itermax = itermax
    self.decrease = decrease
  
  def __call__(self, indiv):
    from random import shuffle

    # computes original fitness.
    if not hasattr( indiv, "fitness" ):
      indiv.fitness = self.evaluation(indiv)

    # dummy class for comparing fitnesses.
    class Dummy():
      def __init__(self, fitness): self.fitness = fitness
   

    indices = range( len(indiv.genes) )
    if self.decrease: indices = [i for i in indices if indiv.genes[i]] 
    iter = 0
    while self.itermax < 0 or iter < self.itermax:

      if len(indices) < 2: break
      shuffle(indices) # random list of genetic indices.

      moved = False
      for j, i in enumerate(indices): 
        indiv.genes[i] = not indiv.genes[i]
        if self.darwin.taboo(self.darwin, indiv ):
          indiv.genes[i] = not indiv.genes[i]
          continue

        new_fitness = Dummy(self.evaluation(indiv))
        iter += 1
        if self.darwin.cmp_indiv(new_fitness, indiv) <= 0: 
          indiv.fitness = new_fitness.fitness
          moved = True
          indices.pop(j)
        else: indiv.genes[i] = not indiv.genes[i]
        
        if iter >= self.itermax: break

      if not moved: break # local minima

    return indiv

          

