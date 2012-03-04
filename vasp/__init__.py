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
__all__ = ['Vasp', 'Extract', 'Specie']
from .extract import Extract
from .specie import Specie
from .functional import Functional as Vasp

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
  from ..misc import read_input
  from . import specie
  from .functional import Functional
  from .extract import Extract

  # names we need to create input.
  input_dict = {"Vasp": Functional, "U": specie.U, "nlep": specie.nlep, 'Extract': Extract}
  if namespace is not None: input_dict.update(namespace)
  return read_input(filepath, input_dict)

def write_incar(path, functional, structure, comm=None):
  """ Writes incar file. """
  from lada import default_comm
  from ..misc.relativepath import RelativePath

  # check what type path is.
  # if not a file, opens one an does recurrent call.
  if path is None: path = "INCAR"
  if not hasattr(path, "write"):
    with open(RelativePath(path).path, "w") as file:
      write_incar(file, functional, structure, comm=comm)
    return

  if comm is None: comm = default_comm
  for line in functional.incar_lines(structure=structure, vasp=functional, comm=comm):
    path.write(line)

def read_incar(path, functional, structure):
  """ Reads from INCAR file and modifies functional accordingly.
  
      Since the functional is meant to work for more than one structure, we
      need the information from the structure this incar refers to.
  """
  from ..misc.relativepath import RelativePath

  # check what type path is.
  # if not a file, opens one an does recurrent call.
  if path is None: path = "INCAR"
  if not hasattr(path, "read"):
    with open(RelativePath(path), "r") as file:
      read_incar(file, functional, structure)

  map = {}
  for line in file:
    if '#' in line: line = line[:line.find('#')]
    if '=' not in line: continue
    key = line[:line.find('=')].rstrip().lstrip()
    value = line[line.find('=')+1:].rstrip().lstrip()
    if '#' in value: value = value[:value.find('#')]
    map[key] = value

  keys = map.keys()
  result = functional.copy()
  for key in (set(keys) & result.params.keys()):
    result.params[key] = map.pop(key)
  for key in (set(keys) & result.params.keys()):
    result.params[key] = map.pop(key)


  


