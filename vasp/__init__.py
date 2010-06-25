""" Module providing an interface to VASP code.

    The interface is separated into 4 conceptual areas:
      - VASP parameterization: mostly contained within L{incar} subpackage.
      - Launching vasp: a single-shot run is performed with L{launch} subpackage.
      - Extracting data from vasp output: to be found in L{extract} subpackage.
      - Methods: such as k-mesh or energy cutoff convergence, strain relaxation....

    The L{Vasp} class  combines the first three concepts together.
    It allows us to launch vasp and retrieve information from the output. It
    checks for errors and avoids running the same job twice. Hence data
    retrieval and vasp calculations can be performed using the same class and
    script. 

    L{version<_vasp.version>} tells for which version of VASP these bindings
    have been compiled.

    L{call_vasp<_vasp.vasp>} allows for direct calls to VASP, with a
    boost.mpi.Communicator as the only argument. VASP input files are expected
    to be found in the current working directory. For an example, see the code
    of L{launch.Launch._run}.
"""
from launch import Launch
from extract import Extract, ExtractGW
from incar import Incar
from kpoints import Density, Gamma
from specie import Specie
from ..opt import __load_vasp_in_global_namespace__
# somehow, mkls don't necessarily get loaded right... Intel linker does not add
# them to libvasp.so. So we have to go through this shit. Same for libsci.so?
if __load_vasp_in_global_namespace__:
  from DLFCN import RTLD_NOW as _RTLD_NOW, RTLD_GLOBAL as _RTLD_GLOBAL
  from sys import getdlopenflags as _getdlopenflags, setdlopenflags as _setdlopenflags
  flags = _getdlopenflags()
  _setdlopenflags(_RTLD_NOW|_RTLD_GLOBAL)
  from _vasp import version, vasp as call_vasp
  import _vasp
  _setdlopenflags(flags)
else:
  import _vasp
call_vasp = _vasp.vasp
version = _vasp.version
""" Vasp version. """
is_vasp_5 = int(version[0]) == 5
""" True if using vasp 5. """
is_vasp_4 = int(version[0]) == 4
""" True if using vasp 4. """
    
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

  def __call__( self, structure, outdir = None, comm = None, repat = None,\
                overwrite=False, keep_calc=False, **kwargs ):
    """ Performs a vasp calculation 
     
        The structure is whatever is given on input. The results are stored in
        directory outdir. The files in L{files.minimal} are copied there, as
        well as any other file named in repat. Other keyword arguements are
        assigned as attributes to a (deepcopy) copy of self prior to actual
        performing the calculation. This way, input parameters can be change
        upon call, while keeping this functor call stateless.

        The return is an L{extract.Extract} object initialized to outdir.

        If successfull results (see L{extract.Extract.success}) already exist
        in outdir, calculations are not repeated. Instead, an extraction object
        for the stored results are given.

        @note: This functor is stateless as long as self and structure can be
               deepcopied correctly.  

        @raise RuntimeError: when computations do not complete.
        @raise IOError: when outdir exists but is not a directory.
    """ 
    from copy import deepcopy
    from os import getcwd
    from os.path import exists, isdir
    from shutil import rmtree
    from boost.mpi import broadcast

    # make this functor stateless.
    this      = deepcopy(self)
    structure = deepcopy(structure)
    outdir    = deepcopy(outdir) if outdir != None else getcwd()
    repat     = deepcopy(repat)  if repat  != None else []
    norun     = kwargs.pop("norun", False)

    is_root = True if comm == None else comm.rank == 0

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value to be changed before launch. 
    for key in kwargs.keys(): setattr(this, key, kwargs[key])

    # First checks if directory outdir exists (and is a directory).
    exists_outdir = broadcast(comm, exists(outdir) if comm.rank == 0 else None, 0) \
                    if comm != None \
                    else exists(outdir)
    if exists_outdir:
      if not overwrite: # check for success
        extract = Extract(comm = comm, directory = outdir)
        if extract.success: return extract # in which case, returns extraction object.
      elif is_root: rmtree(outdir) # overwrite data. 
      if comm != None: comm.barrier() # makes sure directory is not created by other proc!
    
    # Otherwise, performs calculation by calling base class functor.
    super(Vasp, this).__call__( structure=structure, outdir=outdir,\
                                repat=repat, comm=comm, norun=norun, keep_calc=keep_calc )
    
    # checks if result was successful
    extract = Extract(comm = comm, directory = outdir)
    if not norun:
      assert extract.success, RuntimeError("VASP calculation did not complete in %s.\n" % (outdir))

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
