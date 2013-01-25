""" Methods to read structures from file. """
def poscar(path="POSCAR", types=None):
  """ Tries to read a VASP POSCAR file.
      
       :param path: Path to the POSCAR file. Can also be an object with
         file-like behavior.
       :type path: str or file object
       :param types: Species in the POSCAR.
       :type types: None or sequence of str
      
      :return: `pylada.crystal.Structure` instance.
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
      if '#' in line: line = line[:line.find('#')]
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

  def parse_units(line):
    from quantities import a0, meter, centimeter, millimeter, angstrom, emass, \
                           amu, second, millisecond, microsecond, nanosecond,  \
                           picosecond, femtosecond, elementary_charge, coulomb,\
                           hartree, eV, meV, Ry, joule, cal, erg, hertz,       \
                           megahertz, gigahertz, tera, kelvin, newton, dyne,   \
                           h_bar, UnitQuantity, pascal, megapascal, gigapascal,\
                           bar, atm, milli, mol

    auv = UnitQuantity('auv', a0*Ry/h_bar) # velocity
    units = { 'a0': a0, 'bohr': a0, 'm': meter, 'cm': centimeter,
              'mm': millimeter, 'ang': angstrom, 'me': emass, 'amu': amu, 
              's': second, 'ms': millisecond, 'mus': microsecond, 
              'ns': nanosecond, 'ps': picosecond, 'fs': femtosecond,
              'e': elementary_charge, 'c': coulomb, 'hartree': hartree,
              'ha': hartree, 'mha': 1e-3*hartree, 'ev': eV, 'mev': meV,
              'ry': Ry, 'mry': 1e-3*Ry, 'kj': 1e3*joule, 'mol': mol, 
              'kcal': 1e3*cal, 'j': joule, 'erg': erg, 'hz': hertz,
              'mhz': megahertz, 'ghz': gigahertz, 'thz': tera*hertz,
              'k': kelvin, 'n': newton, 'dyne': dyne, 'auv': auv, 'pa': pascal,
              'mpa': megapascal, 'gpa': gigapascal, 'atm': atm, 'bar': bar,
              'atm': atm, 'mbar': milli*bar }
    line = line.replace('cm-1', '1/cm')
    return eval(line, units)


  input = parse_input(file)
  if 'latticecart' in input:
    data = input['latticecart'].splitlines()
    if len(data) == 4:
      units = parse_units(data[0])
      data = data[1:]
    else: units = 1
    cell = array([l.split() for l in data], dtype='float64')
  elif 'latticeabc' in input:
    raise NotImplementedError('Cannot read lattice in ABC format yet.')
  else: 
    raise InputError('Could not find lattice block in input.')
                      
  # create structure
  result = Structure(cell, scale=units)

  # now look for position block.
  units = None
  if 'positionsfrac' in input:
    posdata, isfrac = input['positionsfrac'].splitlines(), True
  elif 'positionsabs' in input:
    posdata, isfrac = input['positionsabs'].splitlines(), False
    try: units = parse_units(posdata[0])
    except: units = None
    else: posdata = posdata[1:]
  else: raise InputError('Could not find position block in input.')
  # and parse it
  for line in posdata:
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




def icsd_cif(filename): 
  """ Reads lattice from the ICSD \*cif files.

      It will not work in the case of other \*cif. 
      It is likely to produce wrong output if the site occupations are fractional.
      If the occupation is > 0.5 it will treat it as 1 and 
      in the case occupation < 0.5 it will treat it as 0 and 
      it will accept all occupation = 0.5 as 1 and create a mess!
  """
  import re
  from os.path import basename
  from numpy.linalg import norm
  from numpy import array, transpose
  from numpy import pi, sin, cos, sqrt, dot

  lines = open(filename,'r').readlines()

  sym_big = 0
  sym_end = 0
  pos_big = 0 
  pos_end = 0

  for l in lines:
      x = l.split()
      if len(x)>0:
          # CELL
          if x[0] == '_cell_length_a':
              if '(' in x[-1]:
                  index = x[-1].index('(')
              else:
                  index = len(x[-1])
              a = float(x[-1][:index])
              
          if x[0] == '_cell_length_b':
              if '(' in x[-1]:
                  index = x[-1].index('(')
              else:
                  index = len(x[-1])
              b = float(x[-1][:index])

          if x[0] == '_cell_length_c':
              if '(' in x[-1]:
                  index = x[-1].index('(')
              else:
                  index = len(x[-1])
              c = float(x[-1][:index])

          if x[0] == '_cell_angle_alpha':
              if '(' in x[-1]:
                  index = x[-1].index('(')
              else:
                  index = len(x[-1])
              alpha = float(x[-1][:index])
              
          if x[0] == '_cell_angle_beta':
              if '(' in x[-1]:
                  index = x[-1].index('(')
              else:
                  index = len(x[-1])
              beta = float(x[-1][:index])

          if x[0] == '_cell_angle_gamma':
              if '(' in x[-1]:
                  index = x[-1].index('(')
              else:
                  index = len(x[-1])
              gamma = float(x[-1][:index])

      # SYMMETRY OPERATIONS

      if len(x)>0 and x[0] == '_symmetry_equiv_pos_as_xyz':
          sym_big = lines.index(l)

      if len(x)>0 and x[0] == '_atom_type_symbol':
          sym_end = lines.index(l)

      # WYCKOFF POSITIONS
      
      if len(x)>0 and x[0] == '_atom_site_attached_hydrogens':
          pos_big = lines.index(l)

      if len(x)>0 and x[0] == '_atom_site_B_iso_or_equiv':
          pos_big = lines.index(l)

      if len(x)>0 and x[0] == '_atom_site_U_iso_or_equiv':
          pos_big = lines.index(l)

      if len(x)>0 and x[0] == '_atom_site_0_iso_or_equiv':
          pos_big = lines.index(l)

      if pos_end == 0 and l in ['\n', '\r\n'] and lines.index(l) > pos_big:
          pos_end = lines.index(l)

  symm_ops = [ '(' + x.split()[1][1:] + x.split()[2] + x.split()[3][:-1] + ')'\
               for x in lines[sym_big+1:sym_end-1] ]

  symm_ops = [re.sub(r'(\d+)', r'\1.', x) for x in symm_ops]
  
  wyckoff = [ [x.split()[0],[x.split()[4],x.split()[5],x.split()[6]],x.split()[7]]\
              for x in lines[pos_big+1:pos_end] ]
  
  wyckoff = [w for w in wyckoff if int(float(w[-1][:4])+0.5) != 0]

  ############## Setting up a good wyckoff list

  for w in wyckoff:
      pom = 0
      for i in range(len(w[0])):
          try:
              int(w[0][i])
              if pom ==0: pom=i
          except:
              pass

      w[0] = w[0][:pom]

      for i in range(3):
          if '(' in w[1][i]: 
              index = w[1][i].index('(')
          else:
              index = len(w[1][i])
          w[1][i] = float(w[1][i][:index])
          
      del w[-1]
  ##########################################
  symbols = list(set([w[0] for w in wyckoff]))
  positions = [[] for i in range(len(symbols))]

  for w in wyckoff:
      symbol = w[0]
      x,y,z = w[1][0],w[1][1],w[1][2]
      for i in range(len(symm_ops)):
          pom = list(eval(symm_ops[i]))
          for j in range(len(pom)):
              if pom[j] <  0.: pom[j] = pom[j]+1.
              if pom[j] >= 0.999: pom[j] = pom[j]-1.

          if not any(norm(array(u)-array(pom)) < 0.01 for u in positions[symbols.index(symbol)]):
              positions[symbols.index(symbol)].append(pom)

  ################ CELL ####################

  a1 = a*array([1.,0.,0.])
  a2 = b*array([cos(gamma*pi/180.),sin(gamma*pi/180.),0.])
  c1 = c*cos(beta*pi/180.)
  c2 = c/sin(gamma*pi/180.)*(-cos(beta*pi/180.)*cos(gamma*pi/180.) + cos(alpha*pi/180.))
  a3 = array([c1, c2, sqrt(c**2-(c1**2+c2**2))])
  cell = array([a1,a2,a3])

  ##########################################

  from lada.crystal import Lattice
  lattice = Lattice()
  lattice.scale = 1.0
  lattice.name = basename(filename)
  lattice.set_cell = transpose(cell)

  for i in range(len(symbols)):
      for j in range(len(positions[i])):
          lattice.add_site = dot(transpose(cell),positions[i][j]), symbols[i]

  lattice.make_primitive()

  return lattice

