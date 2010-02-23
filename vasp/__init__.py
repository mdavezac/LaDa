""" Module providing an interface to VASP code.

    The interface is separated into 4 conceptual areas:
      - VASP parameterization: mostly contained within incar.
      - Launching vasp: a single-shot run is performed with launch.py 
      - Extracting data from vasp output: to be found in extract.py
      - Methods: such as k-mesh or energy cutoff convergence, strain relaxation....

    The Vasp class (in __init__.py) combines the first three concepts together.
    It allows us to launch vasp and retrieve information from the output. It
    checks for errors and avoids running the same job twice. Hence data
    retrieval and vasp calculations can be performed using the same class and
    script). 

    L{version<_vasp.version>} tells for which version of VASP these bindings
    have been compiled.

    L{call_vasp<_vasp.vasp>} allows for direct calls to VASP, with a
    boost.mpi.Communicator as the only argument. VASP input files are expected
    to be found in the current working directory. For an example, see the code
    of L{launch.Launch._run}.
"""
from launch import Launch
from extract import Extract
from incar import Incar
from kpoints import Density, Gamma
from specie import Specie
from _vasp import version, vasp as call_vasp
    
class Vasp(Launch):
  """ Interface to VASP code.
     
      The strength of this interface is that combines vasp calculations, result
      caching, and data retrieval together. 
      
      A vasp run is parameterized using Incar class defined in incar.py.
      It is launched using the Launch class from launch.py class. 
      The results of a successful run is cached in the self.outdir directory. 
      After being launched an object is returned which can extract output data
      from the files in this directory.

      One of the strength of this class is that since results are cached in the
      self.outdir directory, successful calculations are never runned twice.
      This allows us to use the same scripts for generating and retrieving
      data. 
  """

  def __init__(self, *args, **kwargs):
    """ Initializes vasp class. """
    Launch.__init__(self, *args, **kwargs)

  def __call__(self, structure, outdir, mpicomm = None, repat = [], **kwargs):
    """ Performs a vasp calculation 
     
        The structure is whatever is given on input. The results are stored in
        directory outdir. The files in L{files.minimal} are copied there, as
        well as any other file named in repat. Other keyword arguements are
        assigned as attributes to a (deepcopy) copy of self prior to actual
        performing the calculation. This way, input parameters can be change
        upon call, while keeping this functor call stateless.

        The return is an L{extract.Extract} object initialized to outdir.

        If successfull results (see L{extract.Extract.success} and/or
        L{Success<extract._success.Success>}) already exist in outdir,
        calculations are not repeated. Instead, an extraction object for the
        stored results are given.

        @note: This functor is stateless as long as self and structure can be
               deepcopied correctly.  

        @raise RuntimeError: when computations do not complete.
        @raise IOError: when outdir exists but is not a directory.
    """ 
    from copy import deepcopy
    from os.path import exists, isdir

    # make this functor stateless.
    this = deepcopy(self)
    structure = deepcopy(structure)
    outdir = deepcopy(outdir)
    repat = deepcopy(repat)

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value their expected value before launch. 
    if len(kwargs) != 0: 
      for key in kwargs.keys(): getattr(this, key).value = kwargs[key]

    # First checks if directory outdir exists (and is a directory).
    if exists(outdir):
      if not isdir(outdir): raise IOError, "%s exists but is not a directory.\n" % (outdir)
      # checks if it contains a successful run.
      extract = Extract(mpicomm = mpicomm, directory = outdir)
      if extract.success: return extract # in which case, returns extraction object.
    
    # Otherwise, performs calculation by calling base class functor.
    super(Vasp, this).__call__(structure=structure, outdir=outdir, repat=repat, mpicomm=mpicomm)
    
    # checks if result was successful
    extract = Extract(mpicomm = mpicomm, directory = outdir)
    if not extract.success:
      raise RuntimeError, "VASP calculation did not complete in %s.\n" % (outdir)

    return extract

def return_final(looper, *args, **keywords):
  """ Function to pass through a loop, returning final result only. 

      Many methods in this module should be of the form: 
      >>> extract = None
      >>> for extract in method(vasp):
      >>>   # do something
      This function goes through the loop returning only the final result.
      >>> extract = return_final(method.relaxation, structure, vasp, outdir)
  """ 
  result = None
  for result in looper(*args, **keywords): pass
  return result
