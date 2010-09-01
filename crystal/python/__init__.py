""" Contains basic data type and methods for crystal structure and lattices.

    C++ bindings are located in L{_crystal}. A fair number of enhancements are
    added directly within the python code in __init__.py. In practice all
    public interfaces to C++ bindings should be available directly in the
    L{crystal} namespace.
"""
from _crystal import FreezeAtom, which_site, Site, SymmetryOperator, Lattice, to_cartesian,\
                     get_point_group_symmetries, read_structure, sort_layers, \
                     smith_indices, kAtoms, Atom, kAtom, fold_vector, Structure, FreezeCell,\
                     smith_normal_transform,  get_space_group_symmetries, Neighbors, Neighbor, \
                     read_pifile_structure, LayerDepth, to_fractional, linear_smith_index,\
                     rStructure, rAtom, Sites, Atoms, kAtoms # this line not in __all__

from lada.opt.decorators import broadcast_result, add_setter
import A2BX4
import bravais
import gruber
import point_defects

__all__ = [ 'FreezeAtom', 'which_site', 'Sites', 'SymmetryOperator', 'Lattice', 'to_cartesian',\
            'get_point_group_symmetries', 'read_structure', 'sort_layers', 'Site', \
            'smith_indices', 'Atom', 'kAtom', 'fold_vector', 'Structure', 'FreezeCell',\
            'smith_normal_transform', 'get_space_group_symmetries', 'Neighbors', 'Neighbor', \
            'read_pifile_structure', 'LayerDepth', 'to_fractional', 'linear_smith_index', \
            # Below, only true python stuff
            'deform_kpoints', 'read_poscar', 'specie_list', 'write_poscar',\
            'structure_to_lattice', 'fill_structure', 'point_defects', \
            'A2BX4', 'bravais', 'gruber' ]

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
  from numpy import array, dot, transpose
  from . import Structure, Atom
  # checks input
  if check_lattice == True: raise AssertionError, "Not implemented."
  # if types is not none, converts to a list of strings.
  if types != None:
    if isinstance(types, str): types = [types] # can't see another way of doing this...
    elif not hasattr(types, "__getitem__"): types = [str(types)] # single lone vasp.specie.Specie
    else: types = [str(s) for s in types]
      
  if path == None: path = "POSCAR"
  assert exists(path), IOError("Could not find path %s." % (path))
  if isdir(path):
    assert exists(join(path, "POSCAR")), IOError("Could not find POSCAR in %s." % (path))
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
      if types != None:
        assert set(text_types) in set(types) or set(text_types) == set(types), \
               RuntimeError("Unknown species in poscar: %s not in %s." % (set(text_types), set(types)))
      types = text_types
      line = poscar.readline().split()
    assert types != None, RuntimeError("No atomic species given in POSCAR or input.")
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
    
def specie_list(structure):
  """ Returns minimal list of species in alphabetical order. """
  return sorted(list(set([a.type for a in structure.atoms])))

def write_poscar(structure, file, vasp5=False, substitute=None):
  """ Writes a poscar to file. 
  
      >>> with open("POSCAR", "w") as file: write_poscar(structure, file, vasp5=True)

      Species in structures can be substituted for others (when using vasp5 format).
      Below, aluminum atoms are replaced by cadmium atoms. Other atoms are left unchanged.

      >>> with open("POSCAR", "w") as file:
      >>>   write_poscar(structure, file, vasp5=True, substitute={"Al":"Cd"})

      @param structure: The structure to print out.
      @param file: a stream (open file) to which to write.
      @param vasp5: if true, include species in poscar, vasp-5 style.
        Otherwise, does not print specie types.
      @param substitute: if present, will substitute the atom type in the
        structure. Can be incomplete. Only works with vasp5 = True.
      @type substitute: dict
  """
  from numpy import matrix, dot
  file.write(structure.name + "\n")
  file.write(str(structure.scale)+ "\n")
  for i in range(3): file.write("  %f %f %f\n" % tuple(structure.cell[:,i].flat))
  species = specie_list(structure)
  if vasp5: 
    if substitute != None:
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
  

# Adds setter like propertie for easy input 
def _add_atom(which, container):
  """Adds an atom/site to the structure/lattice. 

     >>> structure.add_atom = (0.25,0,0), "Au", 0, FreezeAtom.none
     >>> lattice.add_site = (0.25,0,0), ("Au", "Pd"), 0, FreezeAtom.none
    
     Only the first two arguments are necessary.
     None can be used as a place-holder for the last two arguments.
       - first argument: sequence of three numbers indicating position.
       - second argument: atomic type (or tuple of types for lattices).
       - third argument: index of the site as existing in L{Structure.lattice}
       - fourth argument: whether to freeze
             positional and/or type degrees of freedom.
  """ 
  def _fun(self, args):
    from numpy import array
    args = [x for x in args]
    assert len(args) > 0, RuntimeError("Nothing to set")

    # some skipping around to make sure we parse argument tuple correctly and
    # initialize sites well if type argument is None, or single string.
    if hasattr(args[0], "__iter__"): # assume first argument is a sequence of 3.
      pos = array([x for x in args[0]], dtype="float64")
    else: 
      assert len(args) == 3, RuntimeError("Not sure how to parse these arguments." % (args))
      pos = array([x for x in args], dtype="float64")
      args = args,
    assert pos.size == 3, RuntimeError("First argument is not a sequence of 3." % (args))
    type = None
    if len(args) > 1: type = args[1]
    if type == None:
      result = which()
      result.pos = pos
    else: result = which(pos, type)

    # now everything should be well behaved.
    if len(args) > 2: 
      if args[2] != None: result.site = args[2]
    if len(args) > 3:
      if args[2] != None: result.freeze = args[3]
    getattr(self, container).append(result)
  return _fun

Structure.add_atom = add_setter( _add_atom(Atom, "atoms") )
rStructure.add_k_vec = add_setter( _add_atom(kAtom, "k_vecs") )
rStructure.add_atom = add_setter( _add_atom(rAtom, "atoms") )
Lattice.add_site = add_setter( _add_atom(Site, "sites") )
def _add_atoms(which):
  """ Adds a list of atoms/sites to structure/lattice. 
  
      The argument is a sequence, each item of which could be used with
      L{Structure.add_atom} (or L{Lattice.add_site} when appropriate).

      >>> structure.add_atoms = ((0,0,0), "A"),\\
      >>>                       ((0.25,0,0), "B"),

      Note that the argument must be a sequence of at least 2.

      >>> structure.add_atoms = ((0,0,0), "A") # wrong!
      >>> structure.add_atoms = ((0,0,0), "A"), # Ok

      Items which are None are ignored.
  """
  def _func(self, args):
    for arg in args:
      if arg != None:   setattr(self, which, arg)
  return _func
Structure.add_atoms = add_setter( _add_atoms("add_atom") )
rStructure.add_atoms = add_setter( _add_atoms("add_atom") )
Lattice.add_sites = add_setter( _add_atoms("add_site") )

def _set_types(self, args):
  """ Sets species types in lattice. 

      Sets the species in the lattice using a n-tuple, each item of which is a
      tuple of atomic symbols. The types in each atomic site is set to the
      corresponding item. 

      >>> lattice.set_type = ("In", "Ga"), ("As", "P")

      If an item is None, then that site is unchanged:

      >>> lattice.set_type = ("In", "Ga"), None 

      or 

      >>> lattice.set_type = ("In", "Ga"), 

      will only change the types of the first site.
      Note:

      >>> lattice.set_type = ("In", "Ga") # Error!

      The line above is most likely an error. It will result in setting the
      first site to "In" and second to "Ga".

      >>> lattice.set_type = "In"

      The line above is also an error. It will result in setting the
      first site to "I" and second to "n".
      Similarly, each item must be at least a 2-tuple:

      >>> lattice.set_type = ("In",), # Ok
      >>> lattice.set_type = ("In"), #  Not Okay.

      Finally, the line below, with an empty 2-tuple will change only the second site.

      >>> lattice.set_type = (,), ("As", "P")
      
  """
  assert len(args) <= len(self.sites), RuntimeError("More tuples of types than sites: %s" % (args))
  for types, site in zip(args, self.sites):
    if types == None: continue
    types = [u for u in types if u != None]
    if len(types) == 0: continue
    site.type.clear()
    for type in types: site.type.append(type)

Lattice.set_types = add_setter(_set_types)


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

Structure.set_cell = add_setter(_set_cell)
rStructure.set_cell = add_setter(_set_cell)
Lattice.set_cell = add_setter(_set_cell)

# Makes it easier to set occupation of lattice sites.
def _set_site_type(self, value): 
  """ Sets the type of a site using python intrinsics. """
  self._type.clear()
  if isinstance(value, str): self._type.append(value)
  else: 
    for u in value: self._type.append(u)
def _get_site_type(self): return self._type
Site._type = Site.type 
Site.type = property(_get_site_type, _set_site_type, doc=Site._type.__doc__)
   

# changes __repr__ behavior
def _print_structure(self):
  result  = "# Structure definition.\n"
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
    if atom.site < 0 and atom.freeze == FreezeAtom.none: result += "\n"; continue
    if atom.site < 0: result += "None, "
    else: result += "%i, " % (atom.site)
    result += " %i\n" % (atom.freeze)
  result += "# End of structure definition.\n"
  return result
_print_structure.__doc__ = Structure.__str__.__doc__
Structure.__repr__ = _print_structure

def _print_lattice(self):
  result  = "# Lattice definition.\n"
  result += "lattice = %s()\n" % (self.__class__.__name__)
  result += "lattice.scale = %e\n" % (self.scale)
  result += "lattice.set_cell = (%e, %e, %e),\\\n"\
            "                   (%e, %e, %e),\\\n"\
            "                   (%e, %e, %e)\n"\
            % tuple([x for x in self.cell.flat])
  result += "lattice.name = \"%s\"\n" % (self.name)
  for site in self.sites:
    result += "lattice.add_site = (%e, %e, %e)"  % tuple( x for x in site.pos )
    if len(site.type) == 0: result += "\n"; continue
    if len(site.type) == 1: result += ", \"%s\", " % (site.type[0])
    else:
      result += ", (\"%s\"" % (site.type[0])
      for type in site.type[1:]:  result += ", \"%s\"" % (type)
      result += "), "
    if site.site < 0 and site.freeze == FreezeAtom.none: result += "\n"; continue
    if site.site < 0: result += "None, "
    else: result += "%i, " % (site.site)
    result += " %i\n" % (site.freeze)
  result += "# End of lattice definition.\n"
  return result
_print_lattice.__doc__ = Lattice.__str__.__doc__
Lattice.__repr__ = _print_lattice

def structure_to_lattice(structure):
  """ Converts a structure object to a lattice object. """
  result = Lattice()
  result.cell = structure.cell.copy()
  result.name = structure.name
  result.scale = structure.scale
  for atom in  structure.atoms:
    result.add_site = atom.pos, atom.type
    result.sites[-1].freeze = atom.freeze

  result.make_primitive()
  result.find_space_group()
  return result


def fill_structure(cell, lattice = None):
  """ Returns a structure from knowledge of cell and lattice.

      @param cell: Structure or cell to use to create a complete structure with all atoms.
      @type cell: L{Structure}, L{rStructure}, or numpy 3x3 float64 array
      @param lattice: Back-bone lattice of the super-structure to build. If
        None, will use the *global* lattice set by L{Lattice.set_as_crystal_lattice}.
      @type lattice: L{Lattice}
      @raise RuntimeError: If the filled structure could not be created.
  """
  from _crystal import _fill_structure_impl
  old_lattice = None
  if lattice == None:
    try: Structure().lattice
    except RuntimeError:
      raise RuntimeError("No lattice given on input of fill_structure" +
                         "and global lattice not set either.")
  else: 
    # keeps track of old lattice.
    try: old_lattice = Structure().lattice
    except RuntimeError: pass
    # sets this lattice as the global lattice.
    lattice.set_as_crystal_lattice()

  # creates filled structure.
  result = _fill_structure_impl(cell)

  # Now resets old lattice.
  if old_lattice != None: old_lattice.set_as_crystal_lattice()
  
  return result

