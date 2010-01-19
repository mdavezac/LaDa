""" A GA subpackage defining standard genetic operator for elemental alloys. """

from lada.ga.bitstring import Individual as BitstringIndividual, \
                              Crossover as BitstringCrossover, \
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
    BitstringCrossover.__init__(self, *args, **kwargs)

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

        sStructure().lattice must be set (using  lada.crystal.Lattice.set_as_crystal_lattice(...)).
        @param cell: are the lattice-vectors of the supercell making up the
          elemental superlattice. The epitaxial direction must be given by the
          first column vector.
        @type cell: 3x3 float64 numpy array.




    """
    from lada.crystal import LayerDepth, sort_layers, sStructure, fill_structure
    
    # creates epitaxial structure
    self.structure = sStructure()
    self.structure.cell = cell
    fill_structure(self.structure)
    sort_structure(self.structure)

    ld = LayerDepth(self.structure.cell)
    for i in range(1, len(self.structure.atoms)): 
      if not ld(self.structure.atoms[i-1].pos, self.structure.atoms[i].pos): 
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

class BandgapEvaluator(object):
  """ An evaluator function for bandgaps. """
  def __init__(self, converter, input = "input.xml", lattice = None):
    """ Initializes the bandgap object. 

        @param converter: is a functor which converts between a bitstring and a
          lada.crystal.sStructure, and vice versa.
        @type  converter: duck-typed to L{Converter}
        @param input: is the XML input file for the vff, pescan, 
          (optionally) lattice parameters.
        @type input: string
        @param lattice: Zinc-Blende lattice. If not given, it must be in the
          input file.
        @type lattice: None or lada.crystal.Lattice
    """
    from lada.vff import LayeredVff
    from lada.pescan import BandGap
    from lada import crystal
    import boost.mpi as mpi

    self.converter = converter 
    """ Conversion functor between structures and bitstrings. """

    self.lattice = lattice
    """ Zinc-blend lattice. """
    if self.lattice == None: self.lattice = crystal.Lattice(input)

    self.vff = LayeredVff() # Creates vff functional
    """ Vff functional """
    self.vff.set_mpi(mpi.world) # sets up mpi group for vff.
    self.vff.fromXML(input) # load vff parameters.
    vff.direction = converter.cell[:,0]

    self.escan = BandGap() # creates bandgap functional
    """ Bandgap functional """
    self.escan.set_mpi(mpi.world) # sets mpi group
    self.escan.fromXML(input) # loads escan parameters.
    self.escan.scale = self.lattice.scale

    self.directory_prefix = "indiv"
    """ Directory prefix """
    self.nbcalc = 0
    """ Number of calculations. """

  def __call__(self, indiv):
    """ Computes bandgap of an individual. 
    
        The epitaxial energy is stored in indiv.epi_energy
        The eigenvalues are stored in indiv.eigenvalues
        The VBM and CBM are stored in indiv.bands
        returns the bandgap.
    """
    from boost import mpi
    from lada import crystal
    from lada.opt.changedir import Changedir
    import os
    import shutil
 
    # moves to new calculation directory
    self.bandgap.directory = self.directory_prefix + "_" + str(self.nbcalc)
    if exists(self.bandgap.directory): shutil.rmtree( self.bandgap.directory)
    os.makedirs(self.bandgap.directory)
 
    # moves to calculation directory
    with Changedir(self.bandgap.directory) as pwd:
      # creates structure from bitstring
      self.vff.structure = self.converter(indiv.genes)
      self.vff.init()
     
      # Computes epitaxial structure
      indiv.epi_energy = self.vff.evaluate()
     
      # Computes bandgap
      self.escan.vff_inputfile = "atom_input." + str( mpi.world.rank )
      self.vff.print_escan_input(self.escan.vff_inputfile)
      self.bandgap(self.vff.structure)
      indiv.bands = self.bandgap.bands
 
    # destroy directory if requested.
    if self.bandgap.destroy_directory: shutil.rmtree(self.bandgap.directory)
    
    return indiv.bands.gap


class DipoleEvaluator(BandgapEvaluator):
  """ Evaluates the oscillator strength.

      On top of those quantities saved by base class BandgapEvaluator,
      this class stores the dipole elements in indiv.dipoles.
  """
  def __init__(self, degeneracy = 1e-3, *args, **kwargs): 
    """ Initializes the dipole element evaluator. 

        Dipole elements are evaluated between the VBM and the CBM.
        Bands making up the VBM (CBM) are identified as degenerate if they are
        within "degeneracy" of each other.
    """
    self.degeneracy = degeneracy
    """ Bands making up the VBM (CBM) are identified as degenerate if they are 
        within "degeneracy" of each other. """
    BandgapEvaluator.__init__(self, *args, **kwargs)


  def __call__(self, indiv):
    """ Computes oscillator strength. """
    import os
    import shutil
    import numpy.linalg
    from lada.pescan import dipole_elements
    from lada.opt.changedir import Changedir

    # keeps track of directory destruction
    dodestroy = self.bandgap.escan.destroy_directory 
    self.bandgap.destroy_directory = False
    
    # calls original evaluation function
    super(DipoleEvaluator, self).__call__(self, indiv)

    # moves to calculation directory
    with Changedir(self.bandgap.directory) as pwd:
      # computes dipole elements.
      indiv.dipoles = dipole_elements(self.bandgap, self.degeneracy)

    # destroy calculation directory if requested
    self.bandgap.escan.destroy_directory = dodestroy
    if dodestroy: shutil.rmtree(self.bandgap.directory)

    # return average dipole element.
    result = 0e0
    for u in indiv.dipoles: result += numpy.linalg.norm(u.p)
    result /= float(len(indiv.dipoles))

    return result
