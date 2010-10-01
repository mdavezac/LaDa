""" Contains evaluators for ESCAN properties """
__docformat__ = "restructuredtext en"
from numpy import array as np_array


__all__ = ['cound_calls', 'Bandgap', 'Dipole', 'EffectiveMass']
def count_calls(method):
  """ Increments calls at each call. """
  def wrapped(*args, **kwargs):
    if not hasattr(args[0], "nbcalc"): args[0].nbcalc = 0
    result = method(*args, **kwargs)
    args[0].nbcalc += 1
    return result
  wrapped.__name__ = method.__name__
  wrapped.__doc__  = method.__doc__
  wrapped.__module__ = method.__module__
  wrapped.__dict__.update(method.__dict__)
  return wrapped

  
class Bandgap(object):
  """ An evaluator function for bandgaps at S{Gamma}. """
  def __init__(self, converter, escan, outdir = None, references = None, \
               keep_only_last = True, **kwargs):
    """ Initializes the bandgap object. 

        :Parameters: 
          converter : `lada.ga.escan.elemental.Converter`
            is a functor which converts between a bitstring and a
            `lada.crystal.Structure`, and vice versa.
          escan : `lada.escan.Escan`
            escan functional.
          outdir : str or None
            Directory where to perform calculations. Defaults to None. 
          referemces
	    Reference energies for band-gap calculations. See
            `lada.escan.bandgap`. Defaults to None.
          keep_only_last : boolean
	    If true, only the last calculation is kept. Limits the space
            required for ESCAN optimizations.
          kwargs
            Other keyword arguments to be passed to the bandgap routine.
    """
    from copy import deepcopy
    from ....opt import RelativeDirectory

    super(Bandgap, self).__init__()

    self.converter = converter 
    """ Conversion functor between structures and bitstrings. """
    self.escan = escan
    """ Escan functional """
    self.nbcalc = 0
    """ Number of calculations. """
    self.references = references 
    """ References to compute bandgap.
    
        Can be None (all-electron method), a 2-tuple of reference values
        (folded-spectrum), or a callable returnin a 2-tuple when given a structure.
    """
    self.keep_only_last = keep_only_last
    """ Whethere to keep only last calculations or all. """
    self.kwargs = deepcopy(kwargs)
    """ Additional arguments to be passed to the bandgap functional. """

    self._outdir = RelativeDirectory(path=outdir)
    """ Location of output directory. """

  @property 
  def outdir(self):
    """ Current workding directory. """
    return self._outdir.path
  @outdir.setter 
  def outdir(self, value): self._outdir.path = value

  def __len__(self):
    """ Returns length of bitstring. """
    return len(self.converter)

  def run( self, indiv, outdir = None, comm = None, **kwargs ):
    """ Computes bandgap of an individual. 
    
        :Parameters:
          indiv 
            Individual to compute. Will be converted to a
            `lada.crystal.Structure` using `converter`.
          outdir
            Output directory. 
          comm : boost,mpi.communicator
            MPI communicator.
          kwargs 
            converter, escan, references, outdir  can be overwridden
            on call. This will not affect calls further down the line. Other
            arguments are passed on to the lada.escan.bandgap` function on call.

        :return: an extractor to the bandgap. 

        The epitaxial energy is stored in indiv.epi_energy
        The eigenvalues are stored in indiv.eigenvalues
        The VBM and CBM are stored in indiv.bands
    """
    from shutil import rmtree
    from os.path import join
    from copy import deepcopy
    from boost.mpi import world
    from ....escan import bandgap

    # takes into account input arguments.
    references = kwargs.pop("references", self.references)
    converter  = kwargs.pop( "converter",  self.converter)
    escan      = kwargs.pop(     "escan",      self.escan)
    if comm == None:       comm       = world
    if outdir == None and self.keep_only)last: 
      outdir = self.outdir 
      if comm.rank == 0 and exists(outdir): rmtree(outdir)
      comm.barrier()
    elif outdir == None:
      outdir = join(self.outdir, str(self.nbcalc))
 
    # creates a crystal structure (phenotype) from the genotype.
    structure = converter(indiv.genes)
    # creates argument dictonary to pass on to calculations.
    dictionary = deepcopy(self.kwargs)
    dictionary.update(kwargs) 
    dictionary["comm"]       = comm 
    dictionary["outdir"]     = outdir
    dictionary["workdir"]    = self.outdir
    dictionary["references"] = self.references(structure) if hasattr(references, "__call__")\
                               else references
    if "overwrite" not in dictionary: dictionary["overwrite"]  = True
    # performs calculation.
    out = bandgap(escan, structure, **dictionary)

    # saves some stuff
    indiv.epi_energy = out.energy
    indiv.stress = out.stress.copy()
    indiv.bandgap = out.bandgap
    indiv.vbm = out.vbm
    indiv.cbm = out.cbm
    
    # returns extractor
    return out

  @count_calls
  def __call__(self, *args, **kwargs):
    """ Computes and returns bandgap. 
    
        see self.run(...) for more details.
    """
    return self.run(*args, **kwargs).bandgap

  def __getstate__(self):
    from pickle import dumps
    d = self.__dict__.copy()
    references = self.references
    try:  dumps(d["references"])
    except TypeError: 
      raise RuntimeError("Cannot pickle references in Bandgap evaluator.")
    return d
  def __setstate__(self, arg):
    from pickle import loads
    self.__dict__.update(arg)

  def __repr__(self): 
    """ Returns representation of evaluator. """
    return   "from {0} import {1}\n"\
             "from {2} import {3}\n"\
             "{4}\n\n"\
             "supercell = {5}\n"\
             "converter = {3}(supercell, lattice)\n"\
             "evaluator = {1}(converter, escan_functional)\n"\
             "evaluator.kwargs            = {6}\n"\
             "evaluator.nbcalc            = {7.nbcalc}\n"\
             "evaluator.references        = {7.references}\n"\
             "evaluator.keep_only_last    = {7.keep_only_last}\n"\
             "evaluator.outdir            = {8}\n"\
             .format( self.__class__.__module__, self.__class__.__name__,
                      self.converter.__class__.__module__, self.converter.__class__.__name__,
                      repr(self.escan).replace("functional", "escan_functional"), 
                      repr(self.converter.structure.cell),
                      repr(self.kwargs), self, repr(self._outdir.unexpanded))


class Dipole(Bandgap):
  """ Evaluates the oscillator strength.

      On top of those quantities saved by base class BandgapEvaluator,
      this class stores the dipole elements in indiv.dipoles.
  """
  def __init__(self, *args, **kwargs): 
    """ Initializes the dipole element evaluator.  """
    super(Dipole, self).__init__(*args, **kwargs)

  @count_calls
  def __call__(self, indiv, *args, **kwargs):
    """ Computes the oscillator strength. """
    out = super(Dipole, self).run(indiv, *args, **kwargs)
    indiv.oscillator_strength, indiv.osc_nbstates = out.oscillator_strength()
    return indiv.oscillator_strength


class EffectiveMass(Bandgap):
  """ Evaluates effective mass. """
  def __init__(self, *args, **kwargs): 
    """ Initializes the dipole element evaluator.  """
    from numpy import array
    super(EffectiveMass, self).__init__(*args, **kwargs)

    # isolates effective mass parameters.
    self.emass_dict = { "direction": array([0,0,1]), "order": 2, "nbpoints": None,
                        "stepsize": 1e-2, "lstsq": None }
    for key in self.emass_dict.keys():
      if key in self.kwargs: self.emass_dict[key] = self.kwargs.pop(key)
    # some sanity checks.
    assert self.emass_dict["order"] >= 2

  @count_calls
  def __call__(self, indiv, comm=None, **kwargs):
    """ Computes electronic effective mass. """
    from copy import deepcopy
    from boost.mpi import world
    from ....escan import derivatives

    if comm == None: comm = world

    # isolates effective mass parameters.
    emass_dict = {}
    for key in self.emass_dict.keys():
      emass_dict[key] = kwargs.pop(key, self.emass_dict[key])

    # now gets bandgap.
    out = super(EffectiveMass, self).run(indiv, comm=comm, **kwargs)
    assert out.success, RuntimeError("error in %s" % (out.directory))

    # then prepares parameters for effective mass. 
    # at this point, electronic effective mass, for no good reasons.
    dictionary = deepcopy(self.kwargs)
    dictionary.update(kwargs) 
    dictionary.update(self.emass_dict) 
    dictionary.update(emass_dict) 
    dictionary["outdir"]     = out.directory + "/emass"
    dictionary["eref"]       = out.cbm
    dictionary["structure"]  = out.structure
    # removes band-gap stuff
    dictionary.pop("references", None)
    dictionary.pop("n", None)
    dictionary.pop("overlap_factor", None)
    # at this point, only compute one state.
    dictionary["nbstates"] = 1

    # computes effective mass.
    results = derivatives.reciprocal(self.escan, comm=comm, **dictionary)

    indiv.emass_vbm =  1e0/results[0][2,0]
    return indiv.emass_vbm
