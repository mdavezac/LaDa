""" Contains basic data type and methods for crystal structure and lattices.
    The basic types are imported in this namespace from (and described in)
    module L{_crystal}.
"""
from _crystal import *
import lada.opt.decorators

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
  

@lada.opt.decorators.broadcast_result
def read_poscar(types=None, path=None, check_lattice=False):
  """ Tries to read a VASP POSCAR file,
      
      @param types: species in the POSCAR.
      @type types: sequence of strings or L{lada.vasp.species} or none
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
    if not hasattr(types, "__getitem__"): # probably a lone vasp.specie.Specie instance.
      assert hasattr(types, "symbol"), ValueError("Not sure what argument types is.")
      types = [types.symbol] # makes a list of strings.
    elif isinstance(types, str): types = [types]
    else:
      new_types = []
      for specie in types:
        if hasattr(specie, "symbol"):  # assumes is a vasp.specie.Specie
          new_types.append(specie.symbol) 
        else: new_types.append(specie) # assumes will work as a string.
      types = new_types
      if len(types) == 0: types = None # nothing in sequence

      
  if path == None: path = "POSCAR"
  assert exists(path), "Could not find path %s." % (path)
  if isdir(path):
    assert exists(join(path, "POSCAR")), "Could not find POSCAR in %s." % (path)
    path = join(path, "POSCAR")
  result = Structure()
  with open(path, "r") as poscar:
    # gets name of structure
    result.name = poscar.readline().strip()
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
      if set(text_types) not in set(types):
        raise IOError, "Unknown species in poscar: %s not in %s." % (set(text_types), set(types))
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
    



