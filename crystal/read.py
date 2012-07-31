""" Methods to read structures from file. """
def poscar(path="POSCAR", types=None):
  """ Tries to read a VASP POSCAR file.
      
       :param path: Path to the POSCAR file. Can also be an object with
         file-like behavior.
       :type path: str or file object
       :param types: Species in the POSCAR.
       :type types: None or sequence of str
      
      :return: `lada.crystal.Structure` instance.
  """ 
  import re
  from os.path import join, exists, isdir
  from copy import deepcopy
  from numpy import array, dot, transpose
  from . import Structure

  # if types is not none, converts to a list of strings.
  if types is not None:
    if isinstance(types, str): types = [types] # can't see another way of doing this...
    elif not hasattr(types, "__getitem__"): types = [str(types)] # single lone vasp.specie.Specie
    else: types = [str(s) for s in types]
      
  if path is None: path = "POSCAR"
  if not hasattr(path, 'read'):
    assert exists(path), IOError("Could not find path %s." % (path))
    if isdir(path):
      assert exists(join(path, "POSCAR")), IOError("Could not find POSCAR in %s." % (path))
      path = join(path, "POSCAR")
  result = Structure()
  poscar = path if hasattr(path, "read") else open(path, 'r')
  
  try:
    # gets name of structure
    result.name = poscar.readline().strip()
    if len(result.name) > 0:
      if result.name[0] == "#": result.name = result.name[1:].strip()
    # reads scale
    result.scale = float(poscar.readline().split()[0])
    # gets cell vectors.
    cell = []
    for i in range(3):
      line = poscar.readline()
      assert len(line.split()) >= 3,\
             RuntimeError("Could not read column vector from poscar: %s." % (line))
      cell.append( [float(f) for f in line.split()[:3]] )
    result.cell = transpose(array(cell))
    # checks for vasp 5 input.
    is_vasp_5 = True
    line = poscar.readline().split()
    for i in line: 
      if not re.match(r"[A-Z][a-z]?", i): 
        is_vasp_5 = False
        break
    if is_vasp_5:
      text_types = deepcopy(line)
      if types is not None:
        assert set(text_types) in set(types) or set(text_types) == set(types), \
               RuntimeError( "Unknown species in poscar: {0} not in {1}."\
                             .format(set(text_types), set(types)) )
      types = text_types
      line = poscar.readline().split()
    assert types is not None, RuntimeError("No atomic species given in POSCAR or input.")
    #  checks/reads for number of each specie
    assert len(types) >= len(line), RuntimeError("Too many atomic species in POSCAR.")
    nb_atoms = [int(u) for u in line]
    # Check whether selective dynamics, cartesian, or direct.
    first_char = poscar.readline().strip().lower()[0]
    selective_dynamics = False
    if first_char == 's': 
      selective_dynamics = True
      first_char = poscar.readline().strip().lower()[0]
    # Checks whether cartesian or direct.
    is_direct = first_char not in ['c', 'k']
    # reads atoms.
    for n, type in zip(nb_atoms, types):
      for i in range(n):
        line = poscar.readline().split()
        pos = array([float(u) for u in line[:3]], dtype="float64")
        if is_direct: pos = dot(result.cell, pos)
        result.add_atom(pos=pos, type=type)
        if selective_dynamics:
          for which, freeze in zip(line[3:], ['x', 'y', 'z']):
            if which.lower()[0] == 't':
              result[-1].freeze = getattr(result[-1], 'freeze', '') + freeze
  finally: poscar.close()
            
  return result
    

def castep(file):
  """ Tries to read a castep structure file. """
  from numpy import array, dot
  from ..error import GrepError, IOError
  from ..misc import RelativePath
  from . import Structure
  if isinstance(file, str): 
    if file.find('\n') == -1:
      with open(RelativePath(file).path) as file: return castep(file)
  
  file = [l for l in file]
  # Look for BLOCK lattice_cart
  for i, line in enumerate(file): 
    line = line.split()
    if line[0].lower() == '%block' and line[1].lower() == 'lattice_cart': break
  if i >= len(file) - 1:
    raise GrepError('Could not find lattice_cart block.')

  cell = array( [ file[i+1].split()[:3], file[i+2].split()[:3],
                  file[i+3].split()[:3] ], dtype='float64') 
  result = Structure(cell)

  # now look 
  for i, line in enumerate(file):
    line = line.split()
    if len(line) < 2: continue
    if line[0].lower() == '%block' and line[1].lower() == 'positions_frac':
      break
  if i >= len(file) - 1:
    raise GrepError('Could not find positions_frac block.')
  for line in file[i+1:]:
    line = line.split()
    if len(line) < 2: 
      raise IOError( 'Wrong file format: line with less '                      \
                     'than two items in positions_frac block.')
    if line[0].lower() == '%endblock' and line[1].lower() == 'positions_frac':
      break
    pos = array(line[1:4], dtype='float64')
    pos = dot(result.cell, pos)
    result.add_atom(pos=pos, type=line[0])
  return result
  

