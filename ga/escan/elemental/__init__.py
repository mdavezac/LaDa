""" A GA subpackage defining standard genetic operator for elemental alloys. """
__docformat__ = "restructuredtext en"
from operators import *
from converter import *
from functional import Darwin
from ...bitstring import Individual as BitstringIndividual

class Individual(BitstringIndividual):
  """ An individual for elemental superlattices. 
      
      Comparison between two individuals expect that a bitstring represents an
      elemental superlattice, with translation and inversion symmetries.
  """
  def __init__(self, nmin=0, nmax=20, step=2, dosym=True): 
    """ Initializes a bitstring individual randomly. """
    from random import randint
    import numpy

    assert nmin % 2 == 0 and nmax % 2 == 0 
    self.size = randint(nmin//2, nmax//2) * 2
    super(Individual, self).__init__(self.size)
    self.genes = numpy.array([ int(randint(0,1)) for i in xrange(self.size) ], dtype="int")
    self.dosym = dosym

  def __eq__(self, a): 
    """ Compares two elemental superlattices. """

    if a is None: return False
    if not hasattr(a, "genes"): return False
    if len(a.genes) != len(self.genes): return False
    if getattr(self, 'dosym', False): return all(a.genes - self.genes == 0) 
    
    if self.dosym is False: return all(a.genes - self.genes == 0)
    
    N = len(a.genes)
    for i in range(N):
      # periodicity
      found = True
      for j in range(N):
        if j % dosym != 0: continue
        if a.genes[ (i+j) % N ] != self.genes[j]:
          found = False
          break
      if found: return True
  
      # inversion + periodicity 
      if dosym != 1: continue
      found = True
      for j in range(N):
        if j % dosym != 0: continue
        if a.genes[ -((i+j) % N) ] != self.genes[j]:
          found = False
          break
      if found: return True
  
    return False
