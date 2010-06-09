""" A GA subpackage defining standard genetic operator for elemental alloys. """
__all__ = [ "evaluator" ]
from lada.ga.bitstring import Individual as BitstringIndividual, \
                              Crossover as BitstringCrossover, \
                              Mutation
import evaluator 
class Individual(BitstringIndividual):
  """ An individual for elemental superlattices. 
      
      Comparison between two individuals expect that a bitstring represents an
      elemental superlattice, with translation and inversion symmetries.
  """
  def __init__(self): 
    """ Initializes a bitstring individual randomly. """
    from random import randint
    import numpy

    super(Individual, self).__init__()
    self.genes = numpy.array([ int(randint(0,1)) for i in xrange(self.size) ], dtype="int")
  
  def __eq__(self, a): 
    """ Compares two elemental superlattices. """

    if a == None: return False
    if len(a.genes) != len(self.genes): return False
    
    N = len(a.genes)
    for i in range(N):
      # periodicity
      found = True
      for j in range(N):
        if a.genes[ (i+j) % N ] != self.genes[j]:
          found = False
          break
      if found: return True

      # inversion + periodicity 
      found = True
      for j in range(N):
        if a.genes[ -((i+j) % N) ] != self.genes[j]:
          found = False
          break
      if found: return True

    return False

class Crossover(BitstringCrossover):
  """ A crossover operation between elemental superlattices. """
  def __init__(self, *args, **kwargs):
    BitstringCrossover.__init__(self, *args, **kwargs)

  def __call__(self, a, b):
    """ Performs bitstring crossover between two elemental superlattices.

        Before performing the operation, the origin and direction of one
        superlattice is randomized.
    """
    from random import uniform

    class Dummy: 
      def __init__(self, genes):
        self.genes = b.genes.copy()

    b = Dummy(b)
    if uniform(0,1) < 0.5: b.genes = b.genes[::-1]
    trans = int( uniform(0,1) * len(b.genes) )
    if trans != 0:  
      c = b.genes.copy()
      b.genes[:trans:] = c[len(b.genes)-trans:]
      b.genes[trans:] = c[:len(b.genes)-trans]
    
    return super(Crossover, self).__call__(a,b)


class Converter(object):
  """ Converts a bitstring into an actual superlattice structure, and vice-versa. """

  def __init__(self, cell, lattice = None):
    """ Initializes a functor for bitstring to crystal structure conversion. 

        Structure().L{lattice<crystal.Structure.lattice>} must be set. 
        @param cell: are the lattice-vectors of the supercell making up the
          elemental superlattice. The epitaxial direction must be given by the
          first column vector.
        @type cell: 3x3 float64 numpy array.
        @param lattice: lattice from which to create supercell.
    """
    from lada.crystal import LayerDepth, sort_layers, Structure, fill_structure
    super(Converter, self).__init__()
    
    if lattice != None:
      oldlattice = None
      try: oldlattice = Structure().lattice
      except RuntimeError: pass
      lattice.set_as_crystal_lattice()
    # creates epitaxial tructure
    self.structure = Structure()
    self.structure.cell = cell
    self.structure = sort_layers( fill_structure(self.structure.cell) )

    ld = LayerDepth(self.structure.cell[:,0])
    for i in range(1, len(self.structure.atoms)): 
      if not ld(self.structure.atoms[i-1].pos, self.structure.atoms[i].pos): 
        raise ValueError, "Input structure is not an elemental superlattice."

    if lattice != None:
      if oldlattice != None: oldlattice.set_as_crystal_lattice()
    
  def __call__(self, object):
    """ Conversion function for structures and bitstring. 

        If object is L{crystal.Structure<lada.crystal.Structure>}-like, then converts
        to a bitstring.  Otherwise, expects a bitstring which is converted to
        an L{crystal.Structure<lada.crystal.Structure>}.
    """
    from numpy import array
    from lada.crystal import Structure

    if hasattr(object, "cell"): # crystal structure
      def which(u): # finds type
        return 0 if u.type == self.structure.lattice.sites[u.site].type[0] else 1
      return array([ which(u) for u in object.atoms ])
    else:  # bitstring.
      if len(object) != len(self.structure.atoms): 
        print object
        print len(self.structure.atoms)
        print self.structure
        raise ValueError, "Bitstring and epitaxial structure are not compatible.\n"
      result = Structure(self.structure)
      for i, atom in enumerate(result.atoms):
        atom.type = self.structure.lattice.sites[atom.site].type[int(object[i])]
      assert result.lattice.scale > 0e0
      result.scale = result.lattice.scale
      return result

class LayeredConverter(object):
  """ Converts a bitstring into an actual superlattice structure, and vice-versa. """

  def __init__(self, cell, lattice = None):
    """ Initializes a functor for bitstring to crystal structure conversion. 

        Structure().L{lattice<crystal.Structure.lattice>} must be set. 
        @param cell: are the lattice-vectors of the supercell making up the
          elemental superlattice. The epitaxial direction must be given by the
          first column vector.
        @type cell: 3x3 float64 numpy array.
        @param lattice: lattice from which to create supercell.
    """
    from lada.crystal import LayerDepth, sort_layers, Structure, fill_structure
    super(LayeredConverter, self).__init__()
    
    if lattice != None:
      oldlattice = None
      try: oldlattice = Structure().lattice
      except RuntimeError: pass
      lattice.set_as_crystal_lattice()
    # creates epitaxial tructure
    self.structure = Structure()
    self.structure.cell = cell
    self.structure = sort_layers( fill_structure(self.structure.cell) )

    ld = LayerDepth(self.structure.cell[:,0])
    for i in range(1, len(self.structure.atoms)): 
      if not ld(self.structure.atoms[i-1].pos, self.structure.atoms[i].pos): 
        raise ValueError, "Input structure is not an elemental superlattice."

    if lattice != None:
      if oldlattice != None: oldlattice.set_as_crystal_lattice()
    
  def __call__(self, object):
    """ Conversion function for structures and bitstring. 

        If object is L{crystal.Structure<lada.crystal.Structure>}-like, then converts
        to a bitstring.  Otherwise, expects a bitstring which is converted to
        an L{crystal.Structure<lada.crystal.Structure>}.
    """
    from numpy import array
    from lada.crystal import Structure, LayerDepth

    def generator(struct, all = True):
      layer_depth = LayerDepth(struct.cell[:,0])
      depth, l = None, 0
      for atom in struct.atoms:
        m = layer_depth(atom.pos)
        if depth == None: depth = m  
        elif abs(depth - m) > 1e-12: depth, l = m, l+1
        elif all == False: continue
        yield l, atom

    if hasattr(object, "cell"): # crystal structure
      def which(u): # finds type
        return 0 if u.type == self.structure.lattice.sites[u.site].type[0] else 1
      return array([ which(u) for i, u in generator(object.atoms, False) ])
    else:  # bitstring.
      if len(object) != len([0 for i, u in generator(self.structure.atoms)])
        print object
        print len(self.structure.atoms), len([0 for i, u in generator(self.structure.atoms)])
        print self.structure
        raise ValueError, "Bitstring and epitaxial structure are not compatible.\n"
      result = Structure(self.structure)
      for i, atom in generator(result.atoms, True):
        atom.type = self.structure.lattice.sites[atom.site].type[int(object[i])]
      assert result.lattice.scale > 0e0
      result.scale = result.lattice.scale
      return result
