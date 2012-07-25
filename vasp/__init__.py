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

def parse_incar(path):
  """ Reads INCAR file and returns mapping (keyword, value). """
  from ..error import ValueError
  from ..misc import RelativePath
  if isinstance(path, str): 
    if path.find('\n') == -1:
      with open(RelativePath(path).path) as file: return parse_incar(file)
    else:
      return parse_incar(path.split('\n').__iter__())
  
  lines = []
  for line in path:
    if line.find('#') != -1: line = line[:line.find('#')]
    dummy = [u.lstrip().rstrip() for u in line.split(';')]
    dummy = [u for u in dummy if len(u) > 0]
    if len(dummy) and len(lines) == 0: lines.append(dummy.pop(-1))
    while len(dummy):
      if lines[-1][-1] == '\\':
        lines[-1] = lines[-1][:-1] + dummy.pop(-1)
      else: lines.append(dummy.pop(-1))

  result = {}
  for line in lines:
    if line.find('=') == -1: continue
    keyword, value = [u.rstrip().lstrip() for u in line.split('=')]
    if len(keyword) == 0: raise ValueError('Found empty keword in INCAR.')
    if keyword in result: raise ValueError('Found duplicate keyword {0} in INCAR'.format(keyword))
    result[keyword] = value
  return result
