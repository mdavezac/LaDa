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
    self.evaluated = False
  
  def __eq__(self, a): 
    from math import fabs
    if a == None: return False
    for i, u in enumerate(a.genes):
      if fabs( self.genes[i] - u ) > 1e-18: return False
    return True


class Crossover:
  """ A crossover operation.
  """

  def __init__(self, rate = 0.5):
    self.rate = rate

  def __call__(self, a, b):
    from copy import deepcopy
    from sys import exit
    
    result = deepcopy(a)
    at = self.rate * len(result.genes) 
    result.genes[at:] = b.genes[len(b.genes)-at:]

    return result

class Mutation:
  """ A crossover operation.
  """

  def __init__(self, rate = 0.10):
    self.rate = rate

  def __call__(self, a):
    from copy import deepcopy
    from random import uniform
    
    result = deepcopy(a)

    for i in xrange(len(result.genes)):
      if uniform(0, 1)  < self.rate:
        if result.genes[i] == 1: result.genes[i] = 0
        else: result.genes[i] = 1
    return result



class Mating: 

  """ Chooses between Crossover and Mutation.
  """
  def __init__(self, rate = 0.8, crossover = Crossover(), mutation = Mutation() ):
    self.rate = rate
    self.crossover = crossover
    self.mutation = mutation

  def __call__(self, darwin):
    from random import uniform

    a = darwin.selection(darwin)

    indiv = None
    if uniform(0,1) < self.rate:
      b = a
      while( b == a ): b = darwin.selection(darwin)
      indiv = self.crossover( darwin.population[a], darwin.population[b])
    else: indiv = self.mutation(darwin.population[a])

    if hasattr(indiv, "fitness"): delattr(indiv, "fitness")
    return indiv

