""" Function object to convert to and from bitstrings and nanowires. """
__docformat__ = "restructuredtext en"
class Converter(object):
  """ Converts to and from bitstrings and nanowires. """

  def __init__(self, lattice, growth=(0,0,1), core_radius=2, core_type=0, types=['Si', 'Ge']):
    """ Initializes a multi-shell wire.

        :Parameters:
          lattice = 
          growth 
            Growth direction of the wire.
          core_radius 
            radius of the core. 
          core_type
            Type of the core.
    """
    from copy import deepcopy

    self.lattice = lattice
    """ Zinc-Blende lattice for which the nanowire is constructed. """
    self.growth = growth
    """ Growth direction of the nanowire. """
    self.core_radius = core_radius
    """ Radius of the nanowire's core. """
    self.core_type = core_type
    """ Type of the core. """
    assert self.core_type in [u for a in lattice.sites for u in a.type],\
           ValueError("Could not find atomic-specie of core in lattice.")
    self.types = deepcopy(types)
    """ List of types which GA will optimize.

        Serves to determine what 0 and 1 correspond to in bitstring.
    """
    self.passivant = 'Hg' 
    """ Type of the passivant. """

  def supercell(self, l_plane):
     """ Constructs supercell. """
     from numpy import array, all, abs, round
     from numpy.linalg import norm
     growth = self.growth / norm(self.growth)
     l_plane = int(round(l_plane))
     if all( abs(growth - [0,0,1]) < 1e-12 ):
       # Lattice vector in the column!
       cell = array([[l_plane, 0, 0],[0, l_plane, 0], [0, 0, 1.0]])
       return self.lattice.to_structure(cell, subs={'A':'Hg', 'B':'Hg'})
     if all( abs(growth - [1,1,0]) < 1e-12 ):
       l_xy = int(round(l_plane*sqrt(2)))/2.0
       # Lattice vector in the column!
       cell = array([[-l_xy, 0, 0.5],[l_xy, 0, 0.5], [0, l_plane, 0]]) 
       return self.lattice.to_structure(cell, subs={'A':'Hg', 'B':'Hg'})
     raise NotImplementedError()

  def all_shells(self, bitstring):
    """ Returns an array with all shells, including core. """
    from numpy import concatenate, zeros
    core = zeros(shape=(self.core_radius,), dtype=getattr(bitstring, 'dtype', 'int64'))
    core[:] = self.types.index(self.core_type)
    return concatenate((core, bitstring))

  def superstructure(self, shells):
    """ Returns an empty superstructure of the correct size.
    
        The atoms in the superstructure are correctly centered.
    """
    from numpy import sqrt, dot, floor
    from numpy.linalg import inv
    # constructs the super-structure which will contain the nanowire
    structure = self.supercell(len(shells) * sqrt(2.)/4.0 * 2 + 3)
    
    # Centers atoms in the supercell of the superstructure.
    invcell = inv(structure.cell) 
    for atom in structure.atoms:
      # computes fractional coordinates and recenters atom to middle of cell.
      fractional  = dot(invcell, atom.pos) + [0.5, 0.5, 0] 
      # Move fractional coordinates to [0, 1] interval, taking into account numerical noise.
      fractional -= floor(fractional + 1e-12)
      assert all(fractional + 1e-12 >= 0e0) and all(fractional + 1e-12 < 1e0)
      # Move from fractional back to cartesian coordinates.
      atom.pos = dot(structure.cell, fractional)

    return structure
    
  def to_wire(self, bitstring):
    """ Returns a structure constructed from a bitstring. """
    from numpy import dot
    # creates an array of shell types, including nanowire core.
    shells = self.all_shells(bitstring)
    # create the superstructure which will contain the nanowire.
    nanowire = self.superstructure(shells)

    # finds the center of the cell (and nanowire).
    center = dot(nanowire.cell, [0.5,0.5,0])

    # now iterates over atoms in each shell.
    atomic_shells = self.layer_iterator(nanowire, center, self.growth)
    for shell_index, atoms in enumerate(atomic_shells):
      # finds the type, taking care whether we are inside or outside the nanowire.
      type = self.types[int(shells[shell_index])] if shell_index < len(shells) else self.passivant
      # sets the atoms in that layer to the correct type.
      for atom in atoms: atom.type = type

    return nanowire

  def to_bitstring(self, nanowire):
    """ Constructs a bitstring from a nanowire. """
    from numpy import array, dot

    result = []
    # finds the center of the cell (and nanowire).
    center = dot(nanowire.cell, [0.5,0.5,0])
    # loop over shells in nanowire.
    for layer in self.layer_iterator(nanowire, center, self.growth):
      # gets first atom in layer. 
      atom = layer.next() 
      # break from loop if we reached passivant.
      if atom.type == self.passivant: break 
      # add to list of types.
      result.append( self.types.index(atom.type) )

    return array(result[self.core_radius:])

  def layer_iterator(self, structure, center, direction):
    """ Yields iterators over concentric nanowire shells. 

        Can be use in a double for-loop to initialize or resolve a nano-wire
        type.
    """ 
    from numpy.linalg import norm
    from numpy import dot, sqrt

    # make sure the vector is normalized.
    direction = direction / norm(direction) 
    # layers is a dictionary, where the keys are the shell indices, and
    # the values a list of atoms in the respective shells.
    layers = {}
    for i, atom in enumerate(structure.atoms):
      # vector from center to atom, eg with cartesian coordinates at center.
      pos = atom.pos - center
      # radial distance in cylindrical coordinates, where center is the
      # origin and direction is the z-direction.
      radial_vector = pos - dot(pos, direction) * direction
      radial_norm = norm(radial_vector) * 4e0/sqrt(2e0) # in ML along 011
      # finally, computes shell_index
      shell_index = int(radial_norm+1e-12)
      # add atomic index to shell_index layer.
      if shell_index not in layers: layers[shell_index] = [i]
      else: layers[shell_index].append(i)

    # now loops over shells we just determined, from zero to last.
    for key in sorted(layers.iterkeys()):
      # and yield an iterator over atoms in that shell.
      def iterator():
        """ Iterator over atoms in a single shell. """
        for i in layers[key]: yield structure.atoms[i]
      yield iterator()
  
  def __call__(self, object):
    """ Converts to and from structures and bitstring. """
    if hasattr(object, 'cell'): return self.to_bitstring(object)
    elif hasattr(object, 'genes'): return self.to_wire(object.genes)
    return self.to_wire(object)
 
  def __repr__(self):
    """ Returns string python code representing this object. """
    return "from {0.__class__.__module__} import {0.__class__.__name__}\n"\
           "converter = {0.__class__.__name__}( lattice, growth={1},\n"\
           "                 core_radius={0.core_radius}, core_type='{0.core_type}',\n"\
           "                 types={2} )"\
           .format(self, repr(self.growth), repr(self.types))
