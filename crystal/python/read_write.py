""" Congregates all read/write python routines. """
__docformat__ = "restructuredtext en"
__all__ = ['read_poscar', 'write_poscar', 'write_oldvff', 'read_oldvff', 'icsd_cif']
from lada.opt.decorators import broadcast_result

@broadcast_result(key=True)
def read_poscar(types=None, path=None):
  """ Tries to read a VASP POSCAR file.
      
       :kwarg types: Species in the POSCAR.
       :type types: None or sequence of str
       :kwarg path: Path to the POSCAR file. Can also be an object with
         file-like behavior.
       :type path: str or file object
       :kwarg comm: `mpi.Communicator` over which to read structure.
      
      :return: `lada.crystal.Structure` instance.
  """ 
  import re
  from os.path import join, exists, isdir
  from copy import deepcopy
  from numpy import array, dot, transpose
  from . import Structure, Atom, specie_list

  # if types is not none, converts to a list of strings.
  if types is not None:
    if isinstance(types, str): types = [types] # can't see another way of doing this...
    elif not hasattr(types, "__getitem__"): types = [str(types)] # single lone vasp.specie.Specie
    else: types = [str(s) for s in types]
      
  if path is None: path = "POSCAR"
  assert exists(path), IOError("Could not find path %s." % (path))
  if isdir(path):
    assert exists(join(path, "POSCAR")), IOError("Could not find POSCAR in %s." % (path))
    path = join(path, "POSCAR")
  result = Structure()
  filecontext = path if hasattr(path, "read") else open(path, 'r')
  with filecontext as poscar:
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
    # Checks whether cartesian or direct.
    is_direct = poscar.readline().strip().lower()[0] == "d" 
    # reads atoms.
    for n, type in zip(nb_atoms, types):
      for i in range(n):
        line = poscar.readline().split()
        pos = array([float(u) for u in line[:3]], dtype="float64")
        if is_direct: pos = dot(result.cell, pos)
        result.atoms.append( Atom(pos, type) )
  return result
    
def write_poscar(structure, file, vasp5=False, substitute=None):
  """ Writes a poscar to file. 

      :Parameters:
        structure : `lada.crystal.Structure`
          The structure to print out.
        file : str or file object
          Object with a ''write'' method. If a string, then a file object is
          opened with that filename. The file is overwritten.
        vasp5 : bool
          If true, include species in poscar, vasp-5 style.  Otherwise, does
          not print specie types.
        substitute : None or dict
          If present, will substitute the atom type in the structure. Can be
          incomplete. Only works with vasp5 = True.
  
      >>> with open("POSCAR", "w") as file: write_poscar(structure, file, vasp5=True)

      Species in structures can be substituted for others (when using vasp5 format).
      Below, aluminum atoms are replaced by cadmium atoms. Other atoms are left unchanged.

      >>> with open("POSCAR", "w") as file:
      >>>   write_poscar(structure, file, vasp5=True, substitute={"Al":"Cd"})
  """
  if isinstance(file, str):
    with open(file, 'w') as fileobj: return write_poscar(structure, fileobj, vasp5, substitute)

  from numpy import matrix, dot
  from . import specie_list

  file.write(structure.name + "\n")
  file.write(str(structure.scale)+ "\n")
  for i in range(3): file.write("  %f %f %f\n" % tuple(structure.cell[:,i].flat))
  species = specie_list(structure)
  if vasp5: 
    if substitute is not None:
      for s in species: file.write(" "+ substitute.pop(s,s) +" ")
    else: 
      for s in species: file.write(" "+s+" ")
    file.write("\n")
  for s in species: 
    file.write(" %i " % (len([0 for atom in structure.atoms if atom.type == s])))
  file.write("\nDirect\n")
  inv_cell = matrix(structure.cell).I
  for s in species: 
    for atom in structure.atoms:
      if atom.type != s: continue
      file.write( "  %f %f %f\n" % tuple(dot(inv_cell, atom.pos).flat))
  
def write_oldvff(structure, file, disable_bonds=False):
  """ Writes a structure in the old vff format.

      :Parameters:
        structure : `lada.crystal.Structure`
          The structure to print out.
        file : str or file object
          Object with a ''write'' method. If a string, then a file object is
          opened with that filename. The file is overwritten.
        disable_bonds : bool
          If True, computes bond-types. Otherwise, uses fake bond-types.

      Output is in bhor.
  """
  from quantities import angstrom
  from ..physics import a0, Z
  from . import Neighbors
  if isinstance(file, str):
    with open(file, 'w') as fileobj: return write_oldvff(structure, fileobj)
  
  # gets bond types.
  if not disable_bonds:
    bondtypes = set()
    for atom in structure.atoms:
      for n in Neighbors(structure, 4, atom.pos):
        bondtypes.add( (atom.type, structure.atoms[n.index].type) )
    types = ''
    for bond in bondtypes: 
      types += str(Z(bond[0])) + '0' + str(Z(bond[1])) + ' '
  else: types = '32014 14014 14032 32032'

  file.write( '>lattice-vectors\n'\
              '  {0[0][0]:18.8f}   {0[1][0]:18.8f}  {0[2][0]:18.8f}\n'\
              '  {0[0][1]:18.8f}   {0[1][1]:18.8f}  {0[2][1]:18.8f}\n'\
              '  {0[0][2]:18.8f}   {0[1][2]:18.8f}  {0[2][2]:18.8f}\n\n\n'\
              '>out-put\n'\
              ' {1} {2} \n'\
              '  {0[0][0]:18.8f}   {0[1][0]:18.8f}  {0[2][0]:18.8f}\n'\
              '  {0[0][1]:18.8f}   {0[1][1]:18.8f}  {0[2][1]:18.8f}\n'\
              '  {0[0][2]:18.8f}   {0[1][2]:18.8f}  {0[2][2]:18.8f}\n\n\n'\
              '>atomic positions\n'
              .format(structure.cell, float((structure.scale * angstrom).rescale(a0)), types) )

  for atom in structure.atoms:
    file.write( '  {0.type}   {0.pos[0]:18.8f} {0.pos[1]:18.8f}  {0.pos[2]:18.8f} '\
                '  {0.pos[0]:18.8f} {0.pos[1]:18.8f}  {0.pos[2]:18.8f}\n'\
                .format(atom) )

@broadcast_result(key=True)
def read_oldvff(path):
  """ Tries to read a VASP POSCAR file.
      
      :kwarg path: Path to the POSCAR file. Can also be an
        object with file-like behavior.
      :type path: str or file object.

      :kwarg comm: `mpi.Communicator` over which to read structure.
      
      :return: `lada.crystal.Structure` instance.
  """ 
  from re import compile, M as multline, search
  from numpy import array
  from quantities  import angstrom as AA
  from ..physics import a0
  from . import Structure

  # case where input is a path
  if not hasattr(path, 'read'): 
    with open(path, 'r') as file: return read_oldvff(file)

  # creates both regex expressions.
  re_cell = compile(r'>\s*lattice\s*-?\s*vectors\s*\n'\
                    r'\s*(\S+)\s+(\S+)\s+(\S+)\s*\n'\
                    r'\s*(\S+)\s+(\S+)\s+(\S+)\s*\n'\
                    r'\s*(\S+)\s+(\S+)\s+(\S+)\s*', multline)
  re_pos = compile(r'>\s*atomic\s*-?\s*positions?\s*\n'\
                   r'((?:\s*[A-Z][a-z]?\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*\n)+)',
                   multline)


  # match regex expression to input file.
  string = path.read()
  cell_match = re_cell.search(string)
  pos_match  = re_pos.search(string)
  assert cell_match is not None and pos_match is not None,\
         IOError('File {0} does not seem to be a vff input file.'.format(path.name))

  # creates structure.
  result = Structure()
  result.cell = array(cell_match.groups(), dtype='float64').reshape(3,3).T
  for line in pos_match.group(1).split('\n'):
    data = line.split()
    if len(data) < 4: break
    result.add_atom = data[1:4], data[0]

  scale_match = search(r'>\s*out\s*-?\s*put\s*\n\s*(\S+)', string, multline)
  if scale_match is not None:
    result.scale = float((array(scale_match.group(1), dtype='float64') * a0).rescale(AA))
  else: result.scale = 1 

  return result


@broadcast_result(key=True)
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

      if pos_end == 0 and l == '\n' and lines.index(l) > pos_big:
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

