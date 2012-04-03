""" Module providing an interface to VASP code.

    The interface is separated into 4 conceptual areas:
      - VASP parameterization: mostly contained within L{incar} subpackage.
      - Launching vasp: a single-shot run is performed with L{launch} subpackage.
      - Extracting data from vasp output: to be found in L{extract} subpackage.
      - Methods: such as k-mesh or energy cutoff convergence, strain relaxation....

    The `Vasp` class  combines the first three concepts together.
    It allows us to launch vasp and retrieve information from the output. It
    checks for errors and avoids running the same job twice. Hence data
    retrieval and vasp calculations can be performed using the same class and
    script. 

    `version` tells for which version of VASP these bindings have been
    compiled.
"""
__docformat__ = "restructuredtext en"
__all__ = [ 'GWVasp', 'Vasp', 'Extract', 'Specie',  \
            'Launch', 'Incar', 'RelaxCellShape' ]
from .launch import Launch
from .extract import Extract
from .incar import Incar
from .specie import Specie
from .functional import Functional as Vasp
from .gwfunc import GWFunctional as GWVasp
from .methods import RelaxCellShape

def _vasplib(vasp_library=None):
  """ Returns vasp library string. """
  from lada import vasp_library as global_vasp_lib
  if vasp_library is None: vasp_library = global_vasp_lib
  if vasp_library is None: raise RuntimeError("No vasp library specified.")
  return vasp_library

def version(vasp_library=None, minversion=0):
  """ Vasp version as tuple (major, medium, minor). """
  from _vasp import version as call_version
  result = call_version(_vasplib(vasp_library))
  assert result[0] >= minversion,\
         RuntimeError( "Requested vasp version >= {0}. This is not the case of {1}.\n"\
                       .format(minversion, repr(_vasplib(vasp_library))) )
  return result

def is_vasp_5(vasp_library=None):
  """ True if using vasp 5. """
  try: return version(vasp_library)[0] == 5
  except: return True
def is_vasp_4(vasp_library=None):
  """ True if using vasp 4. """
  try: return version(vasp_library)[0] == 4
  except: return False

def call_vasp(vasp_library=None, comm=None, minversion=0): 
  """ Calls vasp library.  """
  from lada import lada_with_mpi
  from ..mpi import Communicator
  if not lada_with_mpi: raise RuntimeError("Cannot call vasp without mpi.")

  from _vasp import vasp
  comm = Communicator(comm)
  version(vasp_library, minversion)
  vasp(comm, _vasplib(vasp_library))

try: from extract import MassExtract
except ImportError: pass
else: __all__.append(MassExtract.__class__.__name__)
    

def read_input(filepath="input.py", namespace=None):
  """ Specialized read_input function for vasp. 
  
      :Parameters: 
        filepath : str
          A path to the input file.
        namespace : dict
          Addiotional names to include in the local namespace when evaluating
          the input file.

      It add a few names to the input-file's namespace. 
  """
  from ..opt import read_input
  from . import specie
  from .functional import Functional

  # names we need to create input.
  input_dict = { "Vasp": Functional, "GWVasp": GWVasp, "U": specie.U, "nlep": specie.nlep, 
                 "RelaxCellShape": RelaxCellShape, 'Extract': Extract }
  if namespace is not None: input_dict.update(namespace)
  return read_input(filepath, input_dict)
