""" Contains evaluators for Pescan properties """
from numpy import array as np_array
class Bandgap(object):
  """ An evaluator function for bandgaps at S{Gamma}. """
  def __init__(self, converter, input = "input.xml", lattice = None, mpi = None):
    """ Initializes the bandgap object. 

        @param converter: is a functor which converts between a bitstring and a
          lada.crystal.Structure, and vice versa.
        @type  converter: duck-typed to L{Converter}
        @param input: is the XML input file for the vff, pescan, 
          (optionally) lattice parameters.
        @type input: string
        @param lattice: Zinc-Blende lattice. If not given, it must be in the
          input file.
        @type lattice: None or lada.crystal.Lattice
    """
    from boost.mpi import world
    from ....vff import LayeredVff
    from ....escan import BandGap as BGFunctional
    from ... import crystal

    self.world = mpi
    """ MPI Communicator. Defaults to boost.mpi.world. """
    if self.world == None: self.world = world

    self.converter = converter 
    """ Conversion functor between structures and bitstrings. """

    self.lattice = lattice
    """ Zinc-blend lattice. """
    if self.lattice == None: self.lattice = crystal.Lattice(input)
    self.lattice.set_as_crystal_lattice()

    self.vff = LayeredVff(input, self.world) # Creates vff functional
    """ Vff functional """
    self.vff.direction = converter.structure.cell[:,0]
    self.vff.structure.scale = self.lattice.scale

    self.bandgap = BGFunctional(input, self.world) # creates bandgap functional
    """ Bandgap functional """
    self.bandgap.scale = converter.structure

    self.directory_prefix = "indiv"
    """ Directory prefix """
    self.nbcalc = 0
    """ Number of calculations. """


  def __len__(self):
    """ Returns length of bitstring. """
    return len(self.converter.structure.atoms)

  def __call__(self, indiv):
    """ Computes bandgap of an individual. 
    
        The epitaxial energy is stored in indiv.epi_energy
        The eigenvalues are stored in indiv.eigenvalues
        The VBM and CBM are stored in indiv.bands
        returns the bandgap.
    """
    from os.path import exists
    from os import makedirs
    from shutil import rmtree
    from ... import crystal
    from ....opt.changedir import Changedir
 
    # moves to new calculation directory
    self.bandgap.directory = self.directory_prefix + "_" + str(self.nbcalc)
    if self.world.do_print: 
      if exists(self.bandgap.directory): rmtree( self.bandgap.directory)
      makedirs(self.bandgap.directory)
    self.world.barrier()
 
    # moves to calculation directory
    with Changedir(self.bandgap.directory) as pwd:
      # Relax structure
      structure = self.vff( self.converter(indiv.genes) )
      indiv.epi_energy = structure.energy
     
      # Computes bandgap
      self.bandgap.vff_inputfile = "atom_input." + str( self.world.rank )
      self.vff.print_escan_input(self.bandgap.vff_inputfile)
      self.bandgap(structure)
      indiv.bands = self.bandgap.bands
 
    # destroy directory if requested.
    if self.bandgap.destroy_directory: rmtree(self.bandgap.directory)
    
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
    dodestroy = self.bandgap.destroy_directory 
    self.bandgap.destroy_directory = False
    
    # calls original evaluation function
    super(DipoleEvaluator, self).__call__(indiv)

    # moves to calculation directory
    with Changedir(self.bandgap.directory) as pwd:
      # computes dipole elements.
      indiv.dipoles = dipole_elements(self.bandgap, self.degeneracy)

    # destroy calculation directory if requested
    self.bandgap.destroy_directory = dodestroy
    if dodestroy: shutil.rmtree(self.bandgap.directory)

    # return average dipole element.
    result = 0e0
    for u in indiv.dipoles: result += numpy.linalg.norm(u.p)
    result /= float(len(indiv.dipoles))

    return result


class Directness(object):
  """ Objective function for quasi-direct bandgaps. """
  X = np_array( [0,0,1], dtype="float64" )
  """ An X point of the reciprocal zinc-blend lattice. """
  G = np_array( [0,0,0], dtype="float64" )
  """ The S{Gamma} point of the reciprocal zinc-blend lattice. """
  L = np_array( [0.5,0.5,0.5], dtype="float64" )
  """ An L point of the reciprocal zinc-blend lattice. """
  W = np_array( [1, 0.5,0], dtype="float64" )
  """ A W point of the reciprocal zinc-blend lattice. """

  def __init__( self, converter, which, input = "input.xml", lattice = None, mpi = None):
    """ Initializes the Directness object. 

        @param converter: is a functor which converts between a bitstring and a
          lada.crystal.Structure, and vice versa.
        @type  converter: duck-typed to L{Converter}
        @param which: Each element in this sequence contains a k-point, the
          name of that k-point, and a callable returning the reference energy
          as a function of the structure. The first kpoint the one we want the
          lowest.
        @type which: sequence of 3-tuples
        @param input: is the XML input file for the vff, pescan, 
          (optionally) lattice parameters.
        @type input: string
        @param lattice: Zinc-Blende lattice. If not given, it must be in the
          input file.
        @type lattice: None or lada.crystal.Lattice
    """
    from copy import deepcopy
    from numpy.linalg import norm
    from ....vff import LayeredVff
    from ....escan import Escan, method

    self.converter = converter 
    """ Conversion functor between structures and bitstrings. """
    self.which = which
    """ Parameters for computation. """
    assert len(self.which) > 1 
    self.directory_prefix = "indiv"
    """ Directory prefix """
    self.nbcalc = 0

    self.lattice = lattice
    """ Zinc-blend lattice. """
    if self.lattice == None: self.lattice = crystal.Lattice(input)
    self.lattice.set_as_crystal_lattice()

    self.world = mpi
    """ MPI Communicator. Defaults to boost.mpi.world. """
    if self.world == None: self.world = world

    self.vff = LayeredVff(input, self.world) # Creates vff functional
    """ Vff functional """
    self.vff.direction = converter.structure.cell[:,0]

    self.escan = Escan(input, self.world) # creates an escan functional
    """ Bandgap functional """
    # only one state computed.
    self.escan.nbstates = 4 
    # only folded spectrum calculations
    self.escan.method = method.folded


  def __call__(self, indiv):
    """ Evaluates directness of the gap """
    from os.path import exists, join
    from os import makedirs
    from shutil import rmtree
    from numpy import dot as np_dot, matrix as np_matrix
    from numpy.linalg import norm
    from ....escan import Bands, method
    from ....crystal import deform_kpoint
    from ....opt.changedir import Changedir

    assert self.escan.method == method.folded
    assert self.escan.nbstates >= 1

    self.nbcalc += 1
    results = []

    # relax structure.
    structure = self.converter(indiv.genes)
    relaxed, stress = self.vff(structure)
    indiv.epi_energy = relaxed.energy
    indiv.epi_stress = stress

    # create and change directory.
    basedir = self.directory_prefix + "_" + str(self.nbcalc)
    self.world.barrier()
    if self.world.do_print:
      if exists(basedir): rmtree(basedir)

    for kpoint, name, reference in self.which:
      # create new calc directory
      self.escan.directory = join(basedir, name)
      if self.world.do_print: makedirs(self.escan.directory)
      self.world.barrier()

      # change to calc directory
      with Changedir(self.escan.directory) as pwd:

        # computes ideal folded kpoint.
        self.escan.kpoint = deform_kpoint(kpoint, structure.cell, relaxed.cell)
        # gets reference energy from input function.
        self.escan.reference = reference(relaxed)
        #  computes energy
        oldnbstates = self.escan.nbstates 
        if norm(self.escan.kpoint) < 1e-6: self.escan.nbstates /= 2 
        eigenvalues = self.escan(self.vff, relaxed)
        self.escan.nbstates = oldnbstates

        # Finds eigenvalue closes to reference
        mini = eigenvalues[0]
        for eig in eigenvalues[1:]:
          if abs(mini - self.escan.reference) > abs(eig - self.escan.reference):
            mini = eig
        # makes sure all eigenvalues are either above or below references.
        is_good = True
        for eig in eigenvalues[1:]:
          if    (mini - self.escan.reference > 0 and eig - self.escan.reference < 0) \
             or (mini - self.escan.reference < 0 and eig - self.escan.reference > 0):
               is_good = False
               break
        if not is_good: 
          self.escan.method = method.full_diagonalization
          nbvalence = nb_valence_states(structure)
          self.escan.nbstates = nbval + 4
          if norm(self.escan.kpoint) < 1e-6:
            self.escan.nbstates /= 2 
            nbvalence /= 2
          eigenvalues = self.escan(self.vff, relaxed)
          self.escan.nbstates = oldnbstates
          self.escan.method = method.folded
          mini = eig[nbvalence+1]

        # saves eigenvalue result in individual
        setattr(indiv, name + "_eigs", eigenvalues )
        setattr(indiv, name, mini )
        # stores result.
        results.append(mini)

    # destroy directory if requested.
    if self.escan.destroy_directory: rmtree(self.escan.directory)

    # returns maximum offset from gamma 
    gap = results[0] - results[1]
    for e in  results[2:]:
      gap = max( gap, results[0] - e )
    return gap

  def __len__(self):
    """ Returns length of bitstring. """
    return len(self.converter.structure.atoms)

