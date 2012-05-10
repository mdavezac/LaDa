""" A GA subpackage defining standard genetic operator for bitstrings. """
__docformat__ = "restructuredtext en"
from .operators import Crossover, VariableSizeCrossover, Mutation, \
                       GrowthMutation

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
    if size is None: size = self.size
    self.genes = array([ randint(0,1) for i in xrange(size) ])
  
  def __eq__(self, a): 
    from numpy import all, abs
    if a is None: return False
    if len(a.genes) != len(self.genes): return False
    return all(abs(self.genes - a.genes) < 1e-8)

  def __str__(self): return repr(self.genes)
