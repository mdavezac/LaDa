""" Module providing an interface to VASP code.

    The interface is separated into 4 conceptual areas:
      - VASP parameterization: mostly contained within :py:mod:`incar <lada.vasp.incar>` submodule
      - Launching vasp: a single-shot run is performed with a :py:class:`Vasp` object
      - Extracting data from vasp output: to be found in :py:mod:`extract <lada.vasp.extract>` submodule
      - Methods: One can chain vasp runs together for more complex calculations

    The :py:class:`Vasp <lada.vasp.functional.Vasp>` class  combines the first
    three concepts together.  It allows us to launch vasp and retrieve
    information from the output. It checks for errors and avoids running the
    same job twice. Hence data retrieval and vasp calculations can be performed
    using the same class and script. 

    `version` tells for which version of VASP these bindings have been
    compiled.
"""
__docformat__ = "restructuredtext en"
__all__ = ['Vasp', 'Extract', 'Specie', 'MassExtract', 'relax', 'emass', 'read_input', 'exec_input']
from .extract import Extract, MassExtract
from .specie import Specie
from .functional import Vasp
from . import relax, emass

def read_input(filepath="input.py", namespace=None):
  """ Specialized read_input function for vasp. 
  
      :Parameters: 
        filepath : str
          A path to the input file.
        namespace : dict
          Additional names to include in the local namespace when evaluating
          the input file.

      It add a few names to the input-file's namespace. 
  """
  from ..misc import read_input
  from . import specie
  from relax import Epitaxial, Relax

  # names we need to create input.
  input_dict = {"Vasp": Vasp, "U": specie.U, "nlep": specie.nlep, 'Extract': Extract, 
                'Relax': Relax, 'Epitaxial': Epitaxial }
  if namespace is not None: input_dict.update(namespace)
  return read_input(filepath, input_dict)

def exec_input( script, global_dict=None, local_dict=None,
                paths=None, name=None ):
  """ Specialized exec_input function for vasp. """
  from ..misc import exec_input

  # names we need to create input.
  if global_dict is None: global_dict = {}
  for k in __all__:
    if k != 'read_input' and k != 'exec_input': global_dict[k] = globals()[k]
  return exec_input(script, global_dict, local_dict, paths, name)
