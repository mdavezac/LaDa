""" Defines conversion operator from structure to bitstring. """
__docformat__ = "restructuredtext en"
class Converter(object):
  """ Converts to and from bitstrings and nanowires. """

  def __init__(self, lattice, growth=(0,0,1), types=['Si', 'Ge']):
    """ Initializes a multi-shell wire.

        :Parameters:
          lattice 
            Lattice making up the backbone of the superlattice.
          growth : (0,0,1) or (1, 1, 0)
            Growth direction of the wire.
          types : list of str
            Types of atoms in the superlattice. It is used to determine what
            integer correspond to in a bitstring. 
    """
    from copy import deepcopy

    self.lattice = lattice
    """ Zinc-Blende lattice for which the superlattice is constructed. """
    self.growth = growth
    """ Growth direction of the superlattice. """
    self.types = deepcopy(types)
    """ List of types which GA will optimize.

        Serves to determine what 0 and 1 correspond to in bitstring.
    """

  def superstructure(self, n):
    """ Constructs supercell. """
    from numpy import array, all, abs
    if all( abs(array(self.growth) - [0,0,1]) < 1e-12 ):
      assert n % 2 == 0, ValueError("Must have even number of layers. {0}".format(n))
      if n % 4 == 0:
        cell = array([[1, 0.5, 0], [0, 0.5, 0], [0, 0, n//4]], dtype='float64')
      else: 
        cell = array([[1, 0.5, 0], [0, 0.5, 0.5], [0, 0, n//4 + 0.5]], dtype='float64')
    else:
      raise NotImplementedError("Unknown growth direction {0}.".format(self.growth))
    
    return self.lattice.to_structure(cell)
    
  def to_superlattice(self, bitstring):
    """ Returns a structure constructed from a bitstring. """
    from ....crystal import layer_iterator
    # create the superstructure which will contain the superlattices.
    superlattice = self.superstructure(len(bitstring))

    # now iterates over atoms in each layer.
    for shell_index, atoms in enumerate(layer_iterator(superlattice, self.growth)):
      assert shell_index < len(bitstring), (len(bitstring), superlattice.cell)
      assert bitstring[shell_index] in [0, 1], (len(bitstring), superlattice.cell)

      # finds the type, taking care whether we are inside or outside the superlattice.
      type = self.types[int(bitstring[shell_index])]
      # sets the atoms in that layer to the correct type.
      for atom in atoms: atom.type = type

    return superlattice

  def to_bitstring(self, superlattice):
    """ Constructs a bitstring from a superlattice. """
    from numpy import array
    from ....crystal import layer_iterator

    result = []
    # loop over bitstring in superlattice.
    for layer in layer_iterator(superlattice, self.growth):
      # gets first atom in layer. 
      atom = layer.next() 
      # add to list of types.
      result.append( self.types.index(atom.type) )

    return array(result)

  def __call__(self, object):
    """ Converts to and from structures and bitstring. """
    if hasattr(object, 'cell'): return self.to_bitstring(object)
    elif hasattr(object, 'genes'): return self.to_superlattice(object.genes)
    return self.to_superlattice(object)
 
  def __repr__(self):
    """ Returns string python code representing this object. """
    return "from {0.__class__.__module__} import {0.__class__.__name__}\n"\
           "converter = {0.__class__.__name__}( lattice, growth={1}, types={2} )"\
           .format(self, repr(self.growth), repr(self.types))
