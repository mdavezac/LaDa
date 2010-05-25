""" Contains evaluators for Pescan properties """
from numpy import array as np_array
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
  def __init__(self, converter, escan, outdir = None, references = None, **kwargs):
    """ Initializes the bandgap object. 

        @param converter: is a functor which converts between a bitstring and a
          lada.crystal.Structure, and vice versa.
        @type  converter: duck-typed to L{Converter}
        @param escan: escan functional.
        @param **kwargs: Other keyword arguments to be passed to the bandgap routine.
    """
    from copy import deepcopy

    self.converter = converter 
    """ Conversion functor between structures and bitstrings. """
    self.escan = escan
    """ Escan functional """
    self.outdir = "indiv"  if outdir == None else outdir
    """ Output directory.  """
    self.nbcalc = 0
    """ Number of calculations. """
    self.references = references if references != None else lambda x: None
    """ None, or a callback which takes a structure as argument and returns the references. """
    self.kwargs = deepcopy(kwargs)
    """ Additional arguments to be passed to the bandgap functional. """


  def __len__(self):
    """ Returns length of bitstring. """
    return len(self.converter.structure.atoms)

  def run(self, indiv, comm = None):
    """ Computes bandgap of an individual. 
    
        The epitaxial energy is stored in indiv.epi_energy
        The eigenvalues are stored in indiv.eigenvalues
        The VBM and CBM are stored in indiv.bands
        returns an extractor to the bandgap. 
    """
    from os.path import join
    from boost.mpi import world
    from ....escan import bandgap
 
    # performs calculation.
    structure = self.converter(indiv.genes)
    out = bandgap\
          ( 
            self.escan,\
            structure,\
            outdir = join(self.outdir, str(self.nbcalc)),\
            references=self.references(structure),\
            comm = comm if comm != None else world 
          )

    # saves some stuff
    indiv.epi_energy = out.energy
    indiv.stress = out.stress.copy()
    indiv.bandgap = out.bandgap
    indiv.vbm = out.vbm
    indiv.cbm = out.cbm
    
    # returns extractor
    return out

  @count_calls
  def __call__(self, indiv, comm = None):
    """ Computes and returns bandgap. 
    
        see self.run(...) for more details.
    """
    return self.run(indiv, comm).bandgap

class Dipole(Bandgap):
  """ Evaluates the oscillator strength.

      On top of those quantities saved by base class BandgapEvaluator,
      this class stores the dipole elements in indiv.dipoles.
  """
  def __init__(self, *args, **kwargs): 
    """ Initializes the dipole element evaluator.  """
    super(Dipole, self).__init__(*args, **kwargs)

  @count_calls
  def __call__(self, indiv, comm = None):
    """ Computes the oscillator strength. """
    print "nbcalc: ", self.nbcalc
    out = super(Dipole, self).run(indiv, comm)
    indiv.oscillator_strength, indiv.osc_nbstates = out.oscillator_strength()
    return indiv.oscillator_strength
