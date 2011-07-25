""" Crossover operations for CE. """
__docformat__ = "restructuredtext en"
__all__ = ["Crossover", "Evaluator", "Indiv", "Functional", "Mutation"]
from evaluator import Evaluator
from ..bitstring import SwapMutation

class Individual(object):
  """ Individual of a ga for CE. 

      Consists of a bag of ``on`` clusters.
  """
  def __init__(self, maxsize, mean, stddev=None, alwayson=None, alwaysoff=None):
    """ Creates a CE individual. 

        :Parameters:
          maxsize : int
            Maximum size of an individual.
          mean : int 
            Mean size of an individual
          stddev : int or None
            Standard deviation from the mean. 
            If None, then always 
          alwayson : None or sequence
            Clusters that are always in the individual.
          alwaysoff : None or sequence
            Clusters that are never in the individual.
    """ 
    from random import randint
    super(Individual, self).__init__()
    self.maxsize = maxsize
    """ Maximum size of an individual. """

    if stddev == None: stddev = maxsize // 2
    if alwayson == None: alwayson = set()
    if alwaysoff == None: alwaysoff = set()
    genes = set([randint(0, maxsize-1) for u in xrange(randint(max(0, mean-stddev), min(maxsize-1, mean+stddev)))])
    self.genes = list( (genes | alwayson) - alwaysoff)
    """ On clusters. """

  def __eq__(self, other):
    """ True if two individuals are the same. """
    return set(self.genes) == set(other.genes)
  
class Crossover(object):
  """ Crossover operator for CE.

      Creates a new individual from two parents by mixing together only the "on" genes,
      e.g. by taking candy from two bags and putting the into a third. It has
      the advantage that the number of on genes should not grow too fast.
  """
  def __init__(self, rateA=0.5, rateB=None, alwayson=None, alwaysoff=None):
    """ Initializes the crossover operation.

        :Parameters:
          rateA : float 
            Ratio of genes (vs the size of the individual) to take from the first parent.
          rateB : float or None
            Ratio of genes to take from the second parent. If None, then ``rateB = 1-rateA``
          alwayson : None or sequence
            Clusters that are always in the individual.
          alwaysoff : None or sequence
            Clusters that are never in the individual.
    """
    super(Crossover, self).__init__()
    assert rateA > 0 and rateA < 1.0, ValueError("rateA should be in ]0, 1[")

    self.rateA = rateA
    """ How many genes to take from parent A. """
    self.rateB = rateB if rateB != None else 1 - rateA
    """ How many genes to take from parent B. """
    self.alwayson = alwayson if alwayson != None else set()
    """ Clusters that are always in the individual. """
    self.alwaysoff = alwaysoff if alwaysoff != None else set()
    """ Clusters that are never in the individual. """
 
  def __call__(self, a, b):
    """ Create an offspring from two parents. """
    from copy import deepcopy
    from random import shuffle, random

    a,  = list(a.genes), list(b.genes)
    a = shuffle(list(a))[:min(len(a), int(len(a)*random(self.rateA)+0.5))]
    b = shuffle(list(b))[:min(len(b), int(len(b)*random(self.rateB)+0.5))]

    result = deepcopy(a)
    if hasattr(result, "fitness"): delattr(result, "fitness")
    result.genes = list( (set(a+b) | self.alwayson) - self.alwaysoff )

    return result

  def __repr__(self):
    """ Dumps this object to string. """
    string = "{0.__class__.__name__}({0.rateA}".format(self)
    continuous = True
    if abs(1-self.rateB - self.rateA) > 1e-12: 
      string += ", " + str(self.rateB)
    else: continuous = False
    if len(self.alwayson) > 0:
      string += ", "
      if not continuous: string += "alwayson="
      string += repr(list(self.alwayson))
    else: continuous = False
    if len(self.alwaysoff) > 0:
      string += ", "
      if not continuous: string += "alwaysoff="
      string += repr(list(self.alwaysoff))
    string += ")"
    return string

class Mutation(SwapMutation):
  """ Mutation operator for CE.

      Creates a bitstring mutation. 
  """
  def __init__(self, size, rate=0.1, alwayson=None, alwaysoff=None):
    """ Initializes mutation operation. 

        :Parameters:
          size : int
            Number of clusters in the optimization.
          rate : float
            Mutation rate for each bit in the bitstring.
          alwayson : None or sequence
            Clusters which are always on.
          alwaysoff : None or sequence
            Clusters which are always off.
    """
    super(Mutation, self).__init__(rate)
    self.size = size
    """ Size of the bitstring. """
    self.alwayson = alwayson if alwayson != None else set()
    """ Clusters that are always in the individual. """
    self.alwaysoff = alwaysoff if alwaysoff != None else set()
    """ Clusters that are never in the individual. """
 
  def __call__(self, indiv):
    """ Returns invidual with mutated genes. """
    from numpy import array
    from copy import deepcopy

    indiv = deepcopy(indiv)
    indiv.genes = array([(i in indiv.genes) for i in range(self.size)])
    indiv = super(Mutation, self).__init__(indiv)
    indiv.genes = list((set([i for i, u in enumerate(indiv.genes) if u]) | self.alwayson) - self.alwaysoff)
    return indiv

  def __repr__(self):
    """ Dumps this object to string. """
    string = "{0.__class__.__name__}({0.size}, {0.rate}".format(self)
    continuous = True
    if len(self.alwayson) > 0:
      string += ", "
      if not continuous: string += "alwayson="
      string += repr(list(self.alwayson))
    else: continuous = False
    if len(self.alwaysoff) > 0:
      string += ", "
      if not continuous: string += "alwaysoff="
      string += repr(list(self.alwaysoff))
    string += ")"
    return string
