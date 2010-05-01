""" Contains basic data type and methods for crystal structure and lattices.
    The basic types are imported in this namespace from (and described in)
    module L{_crystal}.
"""
from _crystal import *
from lada.opt.decorators import broadcast_result, add_setter

def deform_kpoint(kpoint, ideal, relaxed):
  """ Deform kpoints from ideal cell to relaxed cell. 

      @param kpoint: The kpoint to deform in cartesian coordinates.
      @type kpoint: numpy array
      @param ideal: The original (real-space) cell, as an ideal supercell of the lattice.
      @type ideal: numpy 2d-array
      @param relaxed: The relaxed (real-space) cell.
      @type relaxed: numpy 2d-array
      @return: the kpoint deformed from the ideal reciprocal cell to the
               relaxed reciprocal cell, in cartesian coordinates.
  """
  from numpy import dot, matrix
  from numpy.linalg import inv
  k = dot(ideal.T, kpoint)
  for i in range(3):
    k[i] = float(k[i]) - float( int(k[i]) )
    if k[i] < 0e0: k[i] += 1e0
    if k[i] > 1e0-1e-6: k[i] = 0e0
  return dot(inv(relaxed.T), k.T)
  

@broadcast_result(key=True)
def read_poscar(types=None, path=None, check_lattice=False):
  """ Tries to read a VASP POSCAR file.
      
      @param types: species in the POSCAR.
      @type types: none, or sequence of objects convertible to str 
      @param path: path to the POSCAR file.
      @type path: string
      @param check_lattice: not implemented.
      @type check_lattice: Boolean
      @return: (L{lada.crystal.Structure}) structure on success.
  """ 
  import re
  from os.path import join, exists, isdir
  from copy import deepcopy
  from numpy import array, dot
  from . import Structure, Atom
  # checks input
  if check_lattice == True: raise AssertionError, "Not implemented."
  # if types is not none, converts to a list of strings.
  if types != None:
    if isinstance(types, str): types = [types] # can't see another way of doing this...
    elif not hasattr(types, "__getitem__"): types = [str(types)] # single lone vasp.specie.Specie
    else: types = [str(s) for s in types]
      
  if path == None: path = "POSCAR"
  assert exists(path), "Could not find path %s." % (path)
  if isdir(path):
    assert exists(join(path, "POSCAR")), "Could not find POSCAR in %s." % (path)
    path = join(path, "POSCAR")
  result = Structure()
  with open(path, "r") as poscar:
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
      assert len(line.split()) >= 3, "Could not read column vector from poscar: %s." % (line)
      cell.append( [float(f) for f in line.split()[:3]] )
    result.cell = array(cell)
    # checks for vasp 5 input.
    is_vasp_5 = True
    line = poscar.readline().split()
    for i in line: 
      if not re.match(r"[A-Z][a-z]?", i): 
        is_vasp_5 = False
        break
    if is_vasp_5:
      text_types = deepcopy(line)
      if types != None:
        assert set(text_types) in set(types) or set(text_types) == set(types), \
               IOError("Unknown species in poscar: %s not in %s." % (set(text_types), set(types)))
      types = text_types
      line = poscar.readline().split()
    assert types != None, "No atomic species given in POSCAR or input."
    #  checks/reads for number of each specie
    assert len(types) >= len(line), "Too many atomic species in POSCAR."
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
    

# Adds setter like propertie for easy input 
def _add_atom(which, container):
  """ Defines a property which adds atoms/sites easily. """
  def _add_atom(self, args):
    """ Setter function for the property. """
    from numpy import array
    args = [x for x in args]
    assert len(args) > 0, RuntimeError("Nothing to set")

    result = which()
    if len(args) == 3 and  (not hasattr(args[0], "__iter__")): # assume numpy array only
      result.pos = array(args, dtype="float64") 
    else: # assume that first second argument is position, second is type, etc...
      result.pos = array([x for x in args[0]], dtype="float64")
      if len(args) > 1:
        if args[1] != None: result.type = args[1]
      if len(args) > 2: 
        if args[2] != None: result.site = args[2]
      if len(args) > 3:
        if args[2] != None: result.freeze = args[3]
    getattr(self, container).append(result)
  return _add_atom

Structure.add_atom = add_setter( _add_atom(Atom, "atoms"), 
                                 """Adds an atom to the structure. 

                                    >>> structure.add_atom = (0.25,0,0), "Au", 0, FreezeAtom.none
                                   
                                    Only the first two arguments are necessary.
                                    None can be used as a place-holder for the last two arguments.
                                      - first argument: sequence of three numbers indicating position.
                                      - second argument: atomic type
                                      - third argument: index of the site as existing in L{Structure.lattice}
                                      - fourth argument: whether to freeze
                                            positional and/or type degrees of freedom.
                                 """ )
rStructure.add_k_vec = add_setter( _add_atom(kAtom, "k_vecs"), 
                                   """Adds a k-vector to the structure. 
  
                                      >>> structure.add_atom = (0.25,0,0), "Au", 0+5e0*1j, FreezeAtom.none
                                     
                                      Only the first two arguments are necessary.
                                      None can be used as a place-holder for the last two arguments.
                                        - first argument: sequence of three numbers indicating position.
                                        - second argument: atomic type
                                        - third argument: index of the site as existing in L{Structure.lattice}
                                        - fourth argument: whether to freeze
                                              positional and/or type degrees of freedom.
                                   """ )
rStructure.add_atom = add_setter( _add_atom(rAtom, "atoms"), 
                                 """Adds an atom to the structure. 

                                    >>> structure.add_atom = (0.25,0,0), 1e0, 0, FreezeAtom.none
                                   
                                    Only the first two arguments are necessary.
                                    None can be used as a place-holder for the last two arguments.
                                      - first argument: sequence of three numbers indicating position.
                                      - second argument: atomic type
                                      - third argument: index of the site as existing in L{Structure.lattice}
                                      - fourth argument: whether to freeze
                                            positional and/or type degrees of freedom.
                                 """ )

Lattice.add_site = add_setter( _add_atom(Site, "sites"), 
                               """Adds a site to the lattice. 

                                  >>> structure.add_atom = (0.25,0,0), ("Au", "Pd"), 0, FreezeAtom.none
                                 
                                  Only the first two arguments are necessary.
                                  None can be used as a place-holder for the last two arguments.
                                    - first argument: sequence of three numbers indicating position.
                                    - second argument: atomic type
                                    - third argument: index of the site as existing in L{Structure.lattice}
                                    - fourth argument: whether to freeze
                                          positional and/or type degrees of freedom.
                               """ )


def _set_cell(self, sequence):
  """ Easy input for cells.
  
      >>> structure.cell = (0,0.5,0.5), (0.5,0,0.5), (0.5,0.5,0)
  """
  from numpy import array
  if not hasattr(sequence, "__len__"): # numpy array?
    sequence = [x for x in sequence]
    return _set_cell(self, sequence)

  a0, a1, a2 = [], [], []
  if len(sequence) == 9: a0, a1, a2 = sequence[:3], sequence[3:6], sequence[6:]
  elif len(sequence) == 3:
    a0, a1, a2 = [x for x in sequence[0]], [x for x in sequence[1]], [x for x in sequence[2]]
  else: raise RuntimeError("Don't know how to deal with this much cheese: %s" % (sequence))
  self.cell = array( [a0, a1, a2] )

Structure.set_cell = add_setter(_set_cell, _set_cell.__doc__)
rStructure.set_cell = add_setter(_set_cell, _set_cell.__doc__)
Lattice.set_cell = add_setter(_set_cell, _set_cell.__doc__)

# changes __str__ behavior
def _print_structure(self):
  result  = "# structure definition.\n"
  result += "structure = %s()\n" % (self.__class__.__name__)
  result += "structure.scale = %e\n" % (self.scale)
  result += "structure.set_cell = (%e, %e, %e),\\\n"\
            "                     (%e, %e, %e),\\\n"\
            "                     (%e, %e, %e)\n"\
            % tuple([x for x in self.cell.flat])
  result += "structure.weight = %e\n" % (self.weight)
  result += "structure.name = \"%s\"\n" % (self.name)
  result += "structure.energy = %s\n" % (self.energy)
  for atom in self.atoms:
    result += "structure.add_atom = (%e, %e, %e), \"%s\", " \
              % (atom.pos[0], atom.pos[1], atom.pos[2], str(atom.type))
    if atom.site < 0: result += "None, "
    else: result += "%i, " % (atom.site)
    result += " %i\n" % (atom.freeze)
  return result
_print_structure.__doc__ = Structure.__str__.__doc__
Structure.__str__ = _print_structure

def _print_lattice(self):
  result  = "# lattice definition.\n"
  result += "lattice = %s()\n" % (self.__class__.__name__)
  result += "lattice.scale = %e\n" % (self.scale)
  result += "lattice.set_cell = (%e, %e, %e),\\\n"\
            "                   (%e, %e, %e),\\\n"\
            "                   (%e, %e, %e)\n"\
            % tuple([x for x in self.cell.flat])
  result += "lattice.name = \"%s\"\n" % (self.name)
  for site in self.sites:
    result += "lattice.add_site = (%e, %e, %e) "  % tuple( x for x in site.pos )
    if len(site.type) == 0: result += "\n"; continue
    if len(site.type) == 1: result += ", \"%s\", " % (site.type[0])
    else:
      result += ", (\"%s\"" % (site.type[0])
      for type in site.type[1:]:  result += ", \"%s\"" % (type)
      result += "), "
    if site.site < 0: result += "None, "
    else: result += "%i, " % (site.site)
    result += " %i\n" % (site.freeze)
  return result
_print_lattice.__doc__ = Lattice.__str__.__doc__
Lattice.__str__ = _print_lattice
Lattice.name = ""
