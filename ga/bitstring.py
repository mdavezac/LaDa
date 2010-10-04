""" A GA subpackage defining standard genetic operator for bitstrings. """

class Individual(object):
  """ An individual for bitstrings.

      The size of the individual is given statically by Individual.size.
      @attention: This class does handle mpi at all. Created instances will
        differ from one process to the next.
  """

  size = 10

  def __init__(self, size = None):
    """ Initializes a bitstring individual randomly. 

      @attention: This class does handle mpi at all. Created instances will
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
    self.rate = rate

  def __call__(self, a):
    from random import uniform
    
    for i in xrange(len(a.genes)):
      if uniform(0, 1)  < self.rate:
        if a.genes[i] == 1: a.genes[i] = 0
        else: a.genes[i] = 1
    if hasattr(a, "fitness"): delattr(a, "fitness")
    return a

  def __repr__(self):
    return "from {0} import {1}\nga_operator = {1}({2})"\
           .format(self.__class__.__module__, self.__class__.__name__, self.rate)

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

          

