""" A GA subpackage defining standard genetic operator for bitstrings. """
__docformat__ = "restructuredtext en"

class Individual(object):
  """ An individual for bitstrings.

      The size of the individual is given statically by Individual.size.
      :attention: This class does handle mpi at all. Created instances will
        differ from one process to the next.
  """

  size = 10

  def __init__(self, size = None):
    """ Initializes a bitstring individual randomly. 

      :attention: This class does handle mpi at all. Created instances will
        differ from one process to the next.
    """
    from random import randint
    from numpy import array

    super(Individual, self).__init__()
    if size == None: size = self.size
    self.genes = array([ randint(0,1) for i in xrange(size) ])
  
  def __eq__(self, a): 
    from math import fabs
    if a == None: return False
    for i, u in enumerate(a.genes):
      if fabs( self.genes[i] - u ) > 1e-18: return False
    return True

  def __str__(self): return repr(self.genes)


class VariableSizeCrossover(object):
  """ A crossover operation. """

  def __init__(self, rate = 0.5, Nmin = -1, Nmax = -1):
    """ Initializes crossover. 

        :Parameter:
          Nmin : int
            Minimum size of  a bitstring. Only used if two bitstrings differ in
            size.
          Nmin : int
            Minimum size of  a bitstring. Only used if two bitstrings differ in
            size.
    """
    self.Nmin = Nmin
    """ Minimum bitstring size. """
    self.Nmax = Nmax
    """ Maximum bitstring size. """

  def __call__(self, a, b):
    """ Create new individual from a and b using bitstring crossover. 
    
        :Parameter:
          a 
            First parent individual.
          b 
            This is the other parent individual.

        :return: An offspring individual.
    """
    from copy import deepcopy
    from random import randint

    if len(b.genes) == 0 or len(a.genes) == 0:
      return deepcopy(a if self.rate <= 0.5 else b)
    
    result = deepcopy(a)
    nmin = 0 if self.nmin == -1 else self.nmin
    nmax = min(len(a.indiv.genes), len(b.indiv.genes))
    if self.nmax != -1 and nmax < self.nmax: nmax = self.nmax
    assert nmin < nmax, RuntimeError("Cannot perform crossover.")
    i = randint(nmin, nmax)
    result.genes = result.genes.__class__(list(a.genes[:i]) + list(b.genes[i:]))
    if hasattr(result, "fitness"): delattr(result, "fitness")

    return result
  
  def __repr__(self):
    return "from {0} import {1}\nga_operator = {1}({2.Nmin}, {2.Nmax})"\
           .format(self.__class__.__module__, self.__class__.__name__, self)

class Crossover(object):
  """ A crossover operation. """

  def __init__(self, rate = 0.5):
    self.rate = rate

  def __call__(self, a, b):
    from random import uniform
    
    assert len(a.genes) == len(b.genes), "Bitstring must be of equal lengths"

    i = int( uniform(0, 1) * len(a.genes) )
    j = int( self.rate * len(a.genes) )
    if j >= len(a.genes): 
      a.genes[i:] = b.genes[i:] 
      a.genes[:j-len(a.genes)] = b.genes[:j-len(a.genes)]
    else: a.genes[i:j] = b.genes[i:j]

    if hasattr(a, "fitness"): delattr(a, "fitness")
    return a

  def __repr__(self):
    return "from {0} import {1}\nga_operator = {1}({2})"\
           .format(self.__class__.__module__, self.__class__.__name__, self.rate)

class Mutation(object):
  """ A mutation operation. """

  def __init__(self, rate = 0.10):
    """ Initializes mutation operation. 

        :Parameters:
          rate : float
            Mutation rate for each bit in the bitstring.
    """
    self.rate = rate
    """ Mutation rate for each bitstring. """

  def __call__(self, indiv):
    """ Returns invidual with mutated genes. """
    from random import uniform
    
    rate = self.rate if self.rate >= 0e0 else -float(self.rate)/float(len(indiv.genes))
    for i in xrange(len(indiv.genes)):
      if uniform(0, 1)  < rate:
        indiv.genes[i] = 0 if indiv.genes[i] == 1 else 1 
    if hasattr(indiv, "fitness"): delattr(indiv, "fitness")
    return indiv

  def __repr__(self):
    return "from {0} import {1}\nga_operator = {1}({2})"\
           .format(self.__class__.__module__, self.__class__.__name__, self.rate)
SwapMutation = Mutation
""" Fixed-size bitstring mutation operation. """


class GrowthMutation(object):
  """ Mutation which inserts a bit in the bitstring. """
  def __init__(self, nmin = -1, nmax = -1):
    """ Initiatizes the bit insertion.

        :Parameters:
          nmax : integer
            Maximum size of a bistring. If nmax = -1, there are no limits to the
            size of a bitstring.
          nmin : integer
            Minimum size of a bistring. If nmin = -1, there are no limits to the
            size of a bitstring.
    """
    self.nmax = nmax
    """ Maximum bistring size. """ 
    self.nmin = nmin
    """ Minimum bistring size. """ 
    assert self.nmin < self.nmax or self.nmin == -1 or self.nmax == -1,\
           ValueError("nmin and nmax are incorrect.")

  def __call__(self, indiv):
    """ Inserts extra bit in bitstring. """
    from random import randint, choice
    from itertools import chain

    nmin =  1 if self.nmin != -1 and len(indiv.genes) <= self.nmin else -len(indiv.genes)
    nmax = -1 if self.nmax != -1 and len(indiv.genes) >= self.nmax else len(indiv.genes)
    if nmin > nmax: return indiv

    if nmin < 0 and nmax > 0: 
          n = choice( chain(xrange(nmin, 0), xrange(1, nmax)) )
    else: n = randint(nmin, nmax)
    if n < 0: # negative => remove bit
      l = list(indiv.genes)
      l.pop(n)
    else:     # positive => add bit
      l = list(indiv.genes)
      l.insert(randint[0,2], n)
    indiv.genes = indiv.genes.__class__(l)
     
    if hasattr(indiv, "fitness"): delattr(indiv, "fitness")
    return indiv

  def __repr__(self):
    return "from {0} import {1}\nga_operator = {1}({2.Nmin, 2.Nmax})"\
           .format(self.__class__.__module__, self.__class__.__name__)
    


class LocalSearch(object):
  """ Performs a local search over a bitstring by flipping random bits. """

  def __init__(self, evaluation, darwin, itermax=3, decrease=False):
    """ Initializes a LocalSearch instance.
          - evaluation is a functor or function taking an individual as its argument. 
          - darwin is a class containing a taboo, a selection, and an cmp_indiv
            procedure, as well as a population.
          - itermax is the maximum number of evaluation. 
    """

    assert True, "Not implemented correctly for mpi."
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
        if self.darwin.taboo(indiv ):
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

          

