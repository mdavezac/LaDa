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
  from quantities import angstrom
  from . import Structure

  # if types is not none, converts to a list of strings.
  if types is not None:
    if isinstance(types, str): types = [types] # can't see another way of doing this...
    elif not hasattr(types, "__iter__"): types = [str(types)] # single lone vasp.specie.Specie
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
    result.scale = float(poscar.readline().split()[0]) * angstrom
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
    for n, specie in zip(nb_atoms, types):
      for i in range(n):
        line = poscar.readline().split()
        pos = array([float(u) for u in line[:3]], dtype="float64")
        if is_direct: pos = dot(result.cell, pos)
        result.add_atom(pos=pos, type=specie)
        if selective_dynamics:
          for which, freeze in zip(line[3:], ['x', 'y', 'z']):
            if which.lower()[0] == 't':
              result[-1].freeze = getattr(result[-1], 'freeze', '') + freeze
  finally: poscar.close()
            
  return result
    

def castep(file):
  """ Tries to read a castep structure file. """
  from numpy import array, dot
  from ..periodic_table import find as find_specie
  from ..error import IOError, NotImplementedError, input as InputError
  from ..misc import RelativePath
  from . import Structure
  if isinstance(file, str): 
    if file.find('\n') == -1:
      with open(RelativePath(file).path, 'r') as file: return castep(file)
    else: file = file.splitlines()
  
  file = [l for l in file]

  def parse_input(input):
    """ Retrieves blocks from CASTEP input file. """
    current_block = None
    result = {}
    for line in file:
      if current_block is not None:
        if line.split()[0].lower() == '%endblock': 
          current_block = None
          continue
        result[current_block] += line
      elif len(line.split()) == 0: continue
      elif len(line.split()[0]) == 0: continue 
      elif line.split()[0].lower() == '%block':
        name = line.split()[1].lower().replace('.', '').replace('_', '')
        if name in result:
          raise InputError('Found two {0} blocks in input.'.format(name))
        result[name] = ""
        current_block = name
      else: 
        name = line.split()[0].lower().replace('.', '').replace('_', '')
        if name[-1] in ['=' or ':']: name = name[:-1]
        if name in result:
          raise InputError('Found two {0} tags in input.'.format(name))
        data = line.split()[1:]
        if len(data) == 0: result[name] = None; continue
        if data[0] in [':', '=']: data = data[1:]
        result[name] = ' '.join(data)
    return result


  input = parse_input(file)
  if 'latticecart' in input:
    data = input['latticecart'].splitlines()
    if len(data) == 4:
      raise NotImplementedError('Units cannot be read from CASTEP file')
    cell = array([l.split() for l in data], dtype='float64')
  elif 'latticeabc' in input:
    raise NotImplementedError('Cannot read lattice in ABC format yet.')
  else: 
    raise InputError('Could not find lattice block in input.')
                      
  # create structure
  result = Structure(cell)

  # now look for position block.
  if 'positionsfrac' in input: posdata, isfrac = input['positionsfrac'], True
  elif 'positionsabs' in input: posdata, isfrac = input['positionsabs'], False
  else: raise InputError('Could not find position block in input.')
  # and parse it
  for line in posdata.splitlines():
    line = line.split()
    if len(line) < 2: 
      raise IOError( 'Wrong file format: line with less '                      \
                     'than two items in positions block.')
    pos = array(line[1:4], dtype='float64')
    if isfrac: pos = dot(result.cell, pos)
    try: dummy = int(line[0])
    except: type = line[0]
    else: type = find_specie(atomic_number=dummy).symbol
    result.add_atom(pos=pos, type=type)
    if len(line) == 5: result[-1].magmom = float(line[4])
  return result
  
def crystal(file='fort.34'):
  """ Reads CRYSTAL's external format. """
  from numpy import array, abs, zeros, any, dot
  from numpy.linalg import inv
  from ..crystal import which_site
  from ..misc import RelativePath
  from ..error import IOError
  from ..periodic_table import find as find_specie
  from . import Structure

  if isinstance(file, str):
    if file.find('\n') == -1:
      with open(RelativePath(file).path, 'r') as file: return crystal(file)
    else: file = file.splitlines().__iter__()
  # read first line
  try: line = file.next()
  except StopIteration: raise IOError('Premature end of stream.')
  else: dimensionality, centering, type = [int(u) for u in line.split()[:3]]
  # read cell
  try: cell = array( [file.next().split()[:3] for i in xrange(3)], 
                     dtype='float64' ).T
  except StopIteration: raise IOError('Premature end of stream.')
  result = Structure( cell=cell, centering=centering,
                      dimensionality=dimensionality, type=type, scale=1e0 )
  # read symmetry operators
  result.spacegroup = []
  try: N = int(file.next())
  except StopIteration: raise IOError('Premature end of stream.')
  for i in xrange(N):
    try: op = array( [file.next().split()[:3] for j in xrange(4)],         
                     dtype='float64' )
    except StopIteration: raise IOError('Premature end of stream.')
    else: op[:3] = op[:3].copy().T
    result.spacegroup.append(op)
  result.spacegroup = array(result.spacegroup)

  # read atoms.
  try: N = int(file.next())
  except StopIteration: raise IOError('Premature end of stream.')
  
  for i in xrange(N):
    try: line = file.next().split()
    except StopIteration: raise IOError('Premature end of stream.')
    else: type, pos = int(line[0]), array(line[1:4], dtype='float64')
    if type < 100: type = find_specie(atomic_number=type).symbol
    result.add_atom(pos=pos, type=type, asymmetric=True)

  # Adds symmetrically equivalent structures.
  identity = zeros((4, 3), dtype='float64')
  for i in xrange(3): identity[i, i] == 1
  symops = [u for u in result.spacegroup if any(abs(u - identity) > 1e-8)]
  invcell = inv(result.cell)
  for atom in [u for u in result]:
    for op in symops:
      pos = dot(op[:3], atom.pos) + op[3]
      if which_site(pos, result, invcell=invcell) == -1:
        result.add_atom(pos=pos, type=atom.type, asymmetric=False)

  return result
