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

    self.genome = [ randint(0,1) for i in xrange(size) ]
    self.evaluated = False


class Crossover:
  """ A crossover operation.
  """

  def __init__(self, rate = 0.5):
    self.rate = rate

  def __call__(a, b):
    from copy import deepcopy
    
    result = deepcopy(a)
    at = rate * len(result.genome) 
    result.genome[at:] = b.genome[len(b.genome)-at:]

    return result

class Mutation:
  """ A crossover operation.
  """

  def __init__(self, rate = 0.10):
    self.rate = rate

  def __call__(a):
    from copy import deepcopy
    from random import uniform
    
    result = deepcopy(a)

    for r in result.genome():
      if uniform(0, 1)  < rate:
        if r: r = 0
        else: r = 1
    return result


class Mating: 

  """ Chooses between Crossover and Mutation.
  """
  def __init__(self, rate = 0.2, crossover = Crossover(), mutation = Mutation() ):
    self.rate = 0.8
    self.crossover = crossover
    self.mutation = mutation

  def __call__(self, _selection, _population):

    a, b = _selection(population=_population, size=2)
    indiv = None
    if uniform(0,1) < self.rate: indiv = self.crossover( _population[a], _population[_b]
    else: indiv = self.mutation(_population[a])

    indiv.evaluated = False
    return indiv

