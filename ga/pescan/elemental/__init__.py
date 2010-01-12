""" A GA subpackage defining standard genetic operator for elemental alloys. """

from lada.ga.bistring import Individual as BitstringIndividual, \
                             Crossover as BitstringCrossover,
                             Mutation
class Individual(BitstringIndividual):
  """ An individual for elemental superlattices. 
      
      Comparison between two individuals expect that a bitstring represents an
      elemental superlattice, with translation and inversion symmetries.
  """
  def __init__(self): BitstringIndividual.__init__(self)
  
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
    BitstringCrossover.__init__(self.*args, self.**kwargs)

  def __call__(self, a, b):
    """ Performs bitstring crossover between two elemental superlattices.

        Before performing the operation, the origin and direction of one
        superlattice is randomized.
    """
    from random import uniform

    b = b.copy()
    if uniform(0,1) < 0.5: b = b[::-1]
    trans = int( uniform(0,1) * len(b) )
    if trans != 0:  
      c = b.copy()
      b[:trans:] = c[len(b)-trans:]
      b[trans:] = c[:len(b)-trans]
    
    return BitstringCrossover.__call__(self, a,b)


class Converter(object):
  """ Converts a bitstring into an actual superlattice structure, and vice-versa. """

  def __init__(self, cell):
    """ Initializes a functor for bitstring to crystal structure conversion. 

        @param cell: are the lattice-vectors of the supercell making up the
          elemental superlattice. The epitaxial direction must be given by the
          first column vector.
        @type cell: is an lada.atat.rMatrix3d or derived.

        Prerequisite: sStructure().lattice must be set (using
        lada.crystal.Lattice.set_as_crystal_lattice(...)).
    """
    from lada.crystal import LayerDepth, sort_layers, sStructure, fill_structure
    
    # creates epitaxial structure
    self.structure = sStructure()
    self.structure.cell = cell
    fill_structure(self.structure)
    sort_structure(self.structure)

    ld = LayerDepth(self.structure.cell)
    for i in range(1, len(self.structure.atoms)): 
      if not ld(self.structure.atoms[i-1].pos, self.structure.atoms[i].pos)): 
        raise ValueError, "Input structure is not an elemental superlattice."

    
  def __call__(self, object):
    """ Conversion function for structures and bitstring. 

        If object is lada.crystal.sStructure-like, then converts to a
        bitstring.  Otherwise, expects a bitstring which is converted to an sStructure.
    """
    from lada.crystal import sStructure, Structure

    if hasattr(object, "cell"): # crystal structure
      def which(u): # finds type
        if u.type == self.structure.lattice.sites[u.site].type[0]: 
          return 0
        return 1
      return [ which(u) for u in self.structure.atoms ]
    else:  # bitstring.
      if len(object) != len(self.structure.atoms): 
        raise ValueError, "Bitstring and epitaxial structure are not compatible.\n"
      result = sStructure(self.structure)
      for i, atom in enumerate(result):
        atom.type = self.structure.lattice.sites[atom.site].type[object[i]]
      return Structure(result)

class Bandgap(object):
  """ An evaluator function for bandgaps. """
  def __init__(self, input = "input.xml", lattice = None):
    """ Initializes the bandgap object. 

        @param converter: is a functor which converts between a bitstring and a
          lada.crystal.sStructure, and vice versa.
        @param type: must be a functor with a cell attribute. 
    """
    from lada.vff import LayeredVff
    from lada.pescan import BandGap
    from lada import crystal, atat
    import boost.mpi as mpi

    self.converter = converter 

    self.lattice = lattice
    if self.lattice = None: self.lattice = crystal.lattice(input)

    self.vff = LayeredVff() # Creates vff functional
    self.vff.set_mpi(mpi.world) # sets up mpi group for vff.
    self.vff.fromXML(input) # load vff parameters.
    vff.direction = atat.rVector3d(converter.cell[0,0], converter.cell[1,0], converter.cell[2,0])

    self.escan = BandGap() # creates bandgap functional
    self.escan.set_mpi(mpi.world) # sets mpi group
    self.escan.fromXML(input) # loads escan parameters.
    self.escan.scale = self.lattice.scale

 def __call__(self, indiv):
   """ Computes bandgap of an individual. 
   
       The epitaxial energy is stored in indiv.epi_energy
       The eigenvalues are stored in indiv.eigenvalues
       The VBM and CBM are stored in indiv.bands
       returns the bandgap.
   """
   from boost import mpi
   from lada import crystal

   # creates structure from bitstring
   self.vff.structure = self.converter(indiv.genes)
   self.vff.init()

   # Computes epitaxial structure
   indiv.epi_energy = self.vff.evaluate()

   # Computes bandgap
   self.escan.vff_inputfile = "atom_input." + str( mpi.world.rank )
   self.vff.print_escan_input(self.escan.vff_inputfile)
   self.bandgap.evaluate(self.vff.structure)
   indiv.bands = self.bandgap.bands
   
   return indiv.bands.gap



   

