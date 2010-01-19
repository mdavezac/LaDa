""" Contains evaluators for Pescan properties """
from numpy import array as np_array
class Bandgap(object):
  """ An evaluator function for bandgaps at S{Gamma}. """
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


  def __len__(sefl):
    """ Returns length of bitstring. """
    return len(self.converter.structure.atoms)

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


class Dipole(Bandgap):
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


class Directness(pescan.elemental.evaluator.Bandgap):
  """ Objective function for quasi-direct bandgaps. """
  X = np_array( [0,0,1], dtype="float64" )
  G = np_array( [0,0,0], dtype="float64" )
  L = np_array( [0.5,0.5,0.5], dtype="float64" )
  W = np_array( [1, 0.5,0], dtype="float64" )

  def __init__(self, which = [(Gamma, "Gamma")], *args, **kwargs): 
    super(Eval, self).__init__(self, args, kwargs)
    self.which = which
    
  def __call__(self, indiv):
    """ Evaluates differences between kpoints. """
    from numpy import dot as np_dot
    from lada.escan import Bands

    self.nbcalc += 1
    results = []

    # relax structure.
    structure = converter(indiv)
    self.vff.structure = crystal.Structure(structure)
    self.vff.init()
    indiv.epi_energy = self.vff.evaluate()
    # computes deformation of reciprocal lattice
    deformation = np_dot(self.vff.structure.cell.I.T, structure.cell.T)
    # vff input file
    self.escan.vff_inputfile = "atom_input." + str( mpi.world.rank )

    # create and change directory.
    basedir = self.directory_prefix + "_" + str(self.nbcalc)
    if exists(basedir): shutil.rmtree(basedir)

    for kpoint, name in self.which:
      # create new calc directory
      self.bandgap.directory = join(basedir, name)
      os.makedirs(self.bandgap.directory)

      
      # change to calc directory
      with Changedir(self.bandgap.directory) as pwd:

        # computes ideal folded kpoint.
        self.escan.kpoint = crystal.fold_vector(kpoint, structure.cell.I.T)
        # computes kpoint in deformed structure.
        self.escan.kpoint = np_dot(deformation, self.escan.kpoint)
        # prints escan input
        self.vff.print_escan_input(self.escan.vff_inputfile)
        #  computes bandgap.
        self.bandgap(self.vff.structure)
        # saves bandgap result in individual
        setattr(indiv, name, Bands(self.bandgap.bands) )
        
        results.append( Bands(self.bandgap.bands) )

    # destroy directory if requested.
    if self.bandgap.destroy_directory: shutil.rmtree(self.bandgap.directory)

    # returns maximum offset from gamma 
    result = results[0].cbm - results[1].cbm
    if len(results) > 2: 
      for bands in  results[2:]:
        result = max( result, result[0].cbm - bands.cbm )
    return result
