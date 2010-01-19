""" Module providing an interface to VASP code.

    The interface is separated into 4 conceptual areas:
      - VASP parameterization: mostly contained within incar.py and incar_params.py
      - Launching vasp: a single-shot run is performed with launch.py 
      - Extracting data from vasp output: to be found in extract.py
      - Methods: such as k-mesh or energy cutoff convergence, strain relaxation....

    The Vasp class (in __init__.py) combines the first three concepts together.
    It allows us to launch vasp and retrieve information from the output. It
    checks for errors and avoids running the same job twice. Hence data
    retrieval and vasp calculations can be performed using the same class and
    script). 
"""
from launch import Launch
from extract import Extract
from incar import Incar
from incar_params import *
from kpoints import Density, Gamma
    
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

  def __init__(self):
    """ Initializes vasp class. """
    Launch.__init__(self)

  def __call__(self, structure, outdir, repat = [], **kwargs):
    from copy import deepcopy
    from os.path import exists, isdir

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value their expected value before launch. 
    # A temporary functor object is created, with these attributes modified,
    # and then called.
    if len(kwargs) != 0: 
      this = deepcopy(self)
      for key in kwargs.keys():
        if not hasattr(this, key): raise AttributeError, "Unknown keyword %s.\n" % (key)
        setattr(key, this, kwargs[key]) 
      this(structure, outdir, repat)
    # if no keyword arguments are present, keep going with normal routine.

    # First checks if directory outdir exists (and is a directory).
    if exists(outdir):
      if not isdir(outdir): raise IOError, "%s exists but is not a directory.\n" % (outdir)
      # checks if it contains a successful run.
      extract = Extract(outdir)
      if extract.successful: return extract # in which case, returns extraction object.
    
    # Otherwise, performs calculation by calling base class functor.
    Launch.__call__(self, structure, outdir, repat)
    
    # checks if result was successful
    extract = Extract(outdir)
    if not extract.successful: raise RuntimeError, "VASP calculation did not complete.\n" % (outdir)

    return extract

if __name__ == "__main__":
  from numpy import array as np_array
  from lada import crystal
  from specie import Specie
  
  vasp = Vasp() 
  vasp.species = [Specie("Al", "~/AtomicPotentials/pseudos/K_s")]
  vasp.fftgrid.value = (10,10,10)
  structure = crystal.sStructure()
  structure.scale = 1e0
  structure.cell = np_array([[2,0.5,0.5],[0.0,0,0.5],[0.0,0.5,0]])
  structure.atoms.append( crystal.StrAtom(np_array([0,0,0.0]), "Al") )

  print structure

  vasp(structure)







