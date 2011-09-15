from ._docstring import __doc__ 
__all__ = ['FreezeAtom', 'Atom']

from _crystal import FreezeAtom

def Atom(*args, **kwargs):
  """ Initialize an atom. 

      :Parameters:
        position : list
          Atomic coordinates.  These quantities are accessible as a keyword or
          as the first three arguments of the constructor. 
        type
          For atoms with set or vector types, this is a list of strings representing atomic species.
          For atom with a single string type, this is a single string. Note
          that type can be set using the 4th (to nth for vector and sets)
          arguments rather than this keyword argument. 
        site : int
          Site index. Can only be attained via a keyword only.
        freeze : int
          Site index. Can only be attained via a keyword only.
        atomtype : str
          Defines the atomic type:
            - str: Atomic type is a single string representing a single atomic
                   specie per site.
            - set: Atomic type is a set of strings representing many possible
                   atomic occupation for this site.
            - list: Atomic type is a list of strings representing many
                    possible atomic occupation for this site.

        There are three kinds of atoms. The scalar kind accept only a single
        atomic specie per atomic site, defined in the type attribute. This is
        usefull for defining structures in 90% of applications. However, it may
        be necessary to define a lattice where an atomic site can be occupied
        by any number of species. A lattice is usefull, for instance, in alloy
        applications such as cluster expansion. In that case, two different
        objects are offered, one where the atomic type is list, and the other
        where it is set. 

        There are several ways of constructing atoms, either directly via
        arguments, or via keywords.

        >>> a = Atom(0.25,0.25,0.25, "Au") 
        >>> b = Atom(0.25,0.25,0.25, "Au", "Pd") 
        >>> c = Atom(0.25,0.25,0.25, type="Au") 
        >>> d = Atom(0.25,0.25,0.25, type=["Au", "Pd"]) 
        >>> e = Atom(position=(0.25,0.25,0.25), type=set(["Au", "Pd"])) 

        `a` above is an atom where occupation must be a single atomic specie.
        The position is given by the first three arguments and the occupation
        by the fourth.

        `b` above is an atom where occupation is a *list* of atomic species.
        The position is given by the first three arguments and the occupation
        by the fourth to nth.

        `c` above is an atom where occupation must be a single atomic specie.
        The position is given by the first three arguments and the occupation
        is given by the keyword argument ``type``.

        `d` above is an atom where occupation is a *list* of atomic species.
        The position is given by the first three arguments and the occupation
        is given by the keyword argument ``type``.

        `e` above is an atom where occupation is a *set* of atomic species,
        e.g. a list of unique items, with no duplicates. Both position and type
        are given through keyword arguments.

        Keyword and argument can mixed and matched. 
  """
  from .. import error
  from _crystal import AtomStr, AtomVec, AtomSet, SpecieSet
  position = [0,0,0]
  if len(args) > 0:
    if hasattr(args[0], "__iter__"): 
      if len(args[0]) == 3:
        position = [args[0][0], args[0][1], args[0][2]]
        args = args[1:]
        assert "position" not in kwargs,\
               error.ArgumentError("Position given both in argument and keyword argument.")
      else: raise error.ArgumentError("Unknown argument type {0}.".format(args[0]))
    elif len(args) > 2:
      position = args[:3]
      args = args[3:]
      assert "position" not in kwargs,\
             error.ArgumentError("Position given both in argument and keyword argument.")
    else: raise error.ArgumentError("Incorrect number of arguments in Atom.")
  elif "position" in kwargs: position = kwargs["position"]

  type, atomic_type = "", "scalar"
  if len(args) == 1:
    type, atomic_type = args[0], "scalar"
    assert "type" not in kwargs,\
           error.ArgumentError("Type given both in argument and keyword argument.")
  elif len(args) > 1:
    type, atomic_type = args, "vector"
    assert "type" not in kwargs,\
           error.ArgumentError("Type given both in argument and keyword argument.")
  elif "type" in kwargs:
    type = kwargs["type"]
    if isinstance(type, str): atomic_type = "scalar"
    elif isinstance(type, set) or isinstance(type, SpecieSet): atomic_type = "set"
    else: atomic_type = "list"

  if kwargs.get("atomic_type", atomic_type) != atomic_type:
    if atomic_type == "scalar": atomic_type = kwargs["atomic_type"]
    else: raise error.ArgumentError("Requested a scalar atom but gave {0} atomic species.".format(len(args)))

  if atomic_type == "scalar": atomic_type = AtomStr
  elif atomic_type == "list": atomic_type = AtomVec
  elif atomic_type == "set": atomic_type = AtomSet
  
  args = list(position) + list(type)
  kwargs.pop("position", None)
  kwargs.pop("type", None)
  return atomic_type(*args, **kwargs)





# __all__ = [ 'FreezeAtom', 'which_site', 'Sites', 'SymmetryOperator', 'Lattice', 'to_cartesian',\
#             'get_point_group_symmetries', 'read_structure', 'sort_layers', 'Site', \
#             'smith_indices', 'Atom', 'kAtom', 'fold_vector', 'Structure', 'FreezeCell',\
#             'smith_normal_transform', 'get_space_group_symmetries', 'Neighbors', 'Neighbor', \
#             'read_pifile_structure', 'LayerDepth', 'to_fractional', 'linear_smith_index', \
#             'nb_valence_states', 'to_voronoi', 'gaussian_projector', 'to_cell', 'to_origin', \
#             'is_on_lattice', 'dnc_iterator',
#             # Below, only true python stuff
#             'specie_list', 'read_poscar', 'write_poscar', 'icsd_cif',\
#             'write_oldvff', 'read_oldvff', 'structure_to_lattice', 'fill_structure', \
#             'A2BX4', 'ABX3', 'bravais', 'gruber', 'vasp_ordered', 'binary', 'lattice_context',
#             'layer_iterator', 'equivalence_iterator']
# __docformat__ = "restructuredtext en"

# from _crystal import FreezeAtom, which_site, Site, SymmetryOperator, Lattice, to_cartesian,\
#                      get_point_group_symmetries, read_structure, sort_layers, \
#                      smith_indices, Atom, kAtom, fold_vector, Structure, FreezeCell,\
#                      smith_normal_transform,  get_space_group_symmetries, Neighbors, Neighbor, \
#                      read_pifile_structure, LayerDepth, to_fractional, linear_smith_index,\
#                      nb_valence_states, to_voronoi, to_cell, to_origin, is_on_lattice, \
#                      dnc_iterator,\
#                      rStructure, rAtom, Sites, StringVector # this line not in __all__

# from lada.opt.decorators import add_setter
# from read_write import read_poscar, write_poscar, write_oldvff, read_oldvff, icsd_cif
# import A2BX4
# import ABX3
# import bravais
# import binary
# import gruber
# try: import defects
# except ImportError: pass # required vasp and jobs packages.
# else: __all__.append('defects')

# from contextlib import contextmanager
#   
# def specie_list(structure):
#   """ Returns minimal list of species in alphabetical order. """
#   return sorted(list(set([a.type for a in structure.atoms])))

# # Adds setter like propertie for easy input 
# def _add_atom(which, container):
#   """Adds an atom/site to the structure/lattice. 

#      >>> structure.add_atom = (0.25,0,0), "Au", 0, FreezeAtom.none
#      >>> lattice.add_site = (0.25,0,0), ("Au", "Pd"), 0, FreezeAtom.none
#     
#      Only the first two arguments are necessary.
#      None can be used as a place-holder for the last two arguments.

#        - first argument: sequence of three numbers indicating position.
#        - second argument: atomic type (or tuple of types for lattices).
#        - third argument: index of the site as existing in `Structure.lattice`.
#        - fourth argument: whether to freeze positional and/or type degrees of
#          freedom.

#   """ 

#   def _fun(self, args):
#     from numpy import array
#     args = [x for x in args]
#     assert len(args) > 0, RuntimeError("Nothing to set")

#     # some skipping around to make sure we parse argument tuple correctly and
#     # initialize sites well if type argument is None, or single string.
#     if hasattr(args[0], "__iter__"): # assume first argument is a sequence of 3.
#       pos = array([x for x in args[0]], dtype="float64")
#     else: 
#       assert len(args) == 3, RuntimeError("Not sure how to parse these arguments." % (args))
#       pos = array([x for x in args], dtype="float64")
#       args = args,
#     assert pos.size == 3, RuntimeError("First argument is not a sequence of 3." % (args))
#     type = None
#     if len(args) > 1: type = args[1]
#     if type == None:
#       result = which()
#       result.pos = pos
#     else: result = which(pos, type)

#     # now everything should be well behaved.
#     if len(args) > 2: 
#       if args[2] != None: result.site = args[2]
#     if len(args) > 3:
#       if args[2] != None: result.freeze = args[3]
#     getattr(self, container).append(result)

#   _fun.__doc__ = _add_atom.__doc__
#   return _fun

# Structure.add_atom = add_setter( _add_atom(Atom, "atoms") )
# rStructure.add_k_vec = add_setter( _add_atom(kAtom, "k_vecs") )
# rStructure.add_atom = add_setter( _add_atom(rAtom, "atoms") )
# Lattice.add_site = add_setter( _add_atom(Site, "sites") )
# def _add_atoms(which):
#   """ Adds a list of atoms/sites to structure/lattice. 
#   
#       The argument is a sequence, each item of which could be used with
#       `Structure.add_atom` (or `Lattice.add_site` when appropriate).

#       >>> structure.add_atoms = ((0,0,0), "A"),\\
#       >>>                       ((0.25,0,0), "B"),

#       Note that the argument must be a sequence of at least 2.

#       >>> structure.add_atoms = ((0,0,0), "A") # wrong!
#       >>> structure.add_atoms = ((0,0,0), "A"), # Ok

#       Items which are None are ignored.
#   """
#   def _func(self, args):
#     for arg in args:
#       if arg != None:   setattr(self, which, arg)
#   _func.__doc__ = _add_atoms.__doc__
#   return _func
# Structure.add_atoms = add_setter( _add_atoms("add_atom") )
# rStructure.add_atoms = add_setter( _add_atoms("add_atom") )
# Lattice.add_sites = add_setter( _add_atoms("add_site") )

# def _set_types(self, args):
#   """ Sets species types in lattice. 

#       Sets the species in the lattice using a n-tuple, each item of which is a
#       tuple of atomic symbols. The types in each atomic site is set to the
#       corresponding item. 

#       >>> lattice.set_type = ("In", "Ga"), ("As", "P")

#       If an item is None, then that site is unchanged:

#       >>> lattice.set_type = ("In", "Ga"), None 

#       or 

#       >>> lattice.set_type = ("In", "Ga"), 

#       will only change the types of the first site.
#       Note:

#       >>> lattice.set_type = ("In", "Ga") # Error!

#       The line above is most likely an error. It will result in setting the
#       first site to "In" and second to "Ga".

#       >>> lattice.set_type = "In"

#       The line above is also an error. It will result in setting the
#       first site to "I" and second to "n".
#       Similarly, each item must be at least a 2-tuple:

#       >>> lattice.set_type = ("In",), # Ok
#       >>> lattice.set_type = ("In"), #  Not Okay.

#       Finally, the line below, with an empty 2-tuple will change only the second site.

#       >>> lattice.set_type = (,), ("As", "P")
#       
#   """
#   assert len(args) <= len(self.sites), RuntimeError("More tuples of types than sites: %s" % (args))
#   for types, site in zip(args, self.sites):
#     if types == None: continue
#     types = [u for u in types if u != None]
#     if len(types) == 0: continue
#     site.type.clear()
#     for type in types: site.type.append(type)

# Lattice.set_types = add_setter(_set_types)


# def _set_cell(self, sequence):
#   """ Easy input for cells.
#   
#       >>> structure.cell = (0,0.5,0.5), (0.5,0,0.5), (0.5,0.5,0)
#   """
#   from numpy import array
#   if not hasattr(sequence, "__len__"): # numpy array?
#     sequence = [x for x in sequence]
#     return _set_cell(self, sequence)

#   a0, a1, a2 = [], [], []
#   if len(sequence) == 9: a0, a1, a2 = sequence[:3], sequence[3:6], sequence[6:]
#   elif len(sequence) == 3:
#     a0, a1, a2 = [x for x in sequence[0]], [x for x in sequence[1]], [x for x in sequence[2]]
#   else: raise RuntimeError("Don't know how to deal with this much cheese: %s" % (sequence))
#   dummy = array( [a0, a1, a2], dtype="float64" )
#   self.cell = dummy

# Structure.set_cell = add_setter(_set_cell)
# rStructure.set_cell = add_setter(_set_cell)
# Lattice.set_cell = add_setter(_set_cell)

# # Makes it easier to set occupation of lattice sites.
# def _set_site_type(self, value): 
#   """ Sets the type of a site using python intrinsics. """
#   self._type.clear()
#   if isinstance(value, str): self._type.append(value)
#   else: 
#     for u in value: self._type.append(u)
# def _get_site_type(self): return self._type
# Site.type = property(_get_site_type, _set_site_type, doc=Site._type.__doc__)
#    
# def _repr_site_type(self):
#   "Represents the type of a lattice site."
#   if len(self) == 1: return repr(self[0])
#   return repr(tuple([u for u in self]))
# StringVector.__repr__ = _repr_site_type

# # changes __repr__ behavior
# def _repr_lattstr(self, atoms="atoms", name="structure", add_atoms="add_atoms"):
#   """ Represents lattices and/or structures. """
#   result  = "# {c}{name} definition.\n".format(c=name[0].upper(), name=name[1:])
#   # prints all but atoms.
#   result += "from {0} import {1}\n{name} = {1}()\n"\
#             .format(self.__class__.__module__, self.__class__.__name__, self, name=name)
#   # force representation of these values.
#   do_represent = ['name', 'scale', 'energy', 'weight', 'magmom']
#   width = max(len(u) for u in do_represent)
#   for key in do_represent:
#     if not hasattr(self, key): continue
#     try:
#       r = repr( getattr(self, key) )
#       assert r[0] != '<' and r[-1] != '>'
#     except: result += "# Could not represent {name}.{0}.".format(key, name=name)
#     else:
#       start = "{name}.{0: <{1}} = ".format(key, width, name=name)
#       r = r.replace('\n', '\\\n{0:{1}}'.format(' ', len(start)))
#       result += "{0}{1}\n".format(start, r)
#   # tries to print rest of dictionary.
#   for key, value in self.__dict__.items(): 
#     if key in do_represent: continue
#     try: 
#       r = repr(value)
#       assert r[0] != '<' and r[-1] != '>'
#     except: result += "# Could not represent {name}.{0}.".format(key, name=name)
#     else: result += "{name}.{0} = {1}\n".format(key, repr(value), name=name)

#   # now adds cell
#   width = max(len("{0:.0f}".format(u)) for u in self.cell.flat) 
#   precision = ["{0}".format(u) for u in self.cell.flat]
#   precision = max(len(u[(u.index('.') + 1 if '.' in u else 0):]) for u in precision)
#   width = precision + width
#   result += "{start}({0[0][0]:{1}.{2}f}, {0[0][1]:{1}.{2}f}, {0[0][2]:{1}.{2}f}),\\\n"\
#             "{empty}({0[1][0]:{1}.{2}f}, {0[1][1]:{1}.{2}f}, {0[1][2]:{1}.{2}f}),\\\n"\
#             "{empty}({0[2][0]:{1}.{2}f}, {0[2][1]:{1}.{2}f}, {0[2][2]:{1}.{2}f})\n"\
#             .format(self.cell, width, precision, 
#                     start = "{0}.set_cell = ".format(name),\
#                     empty = "".join(" " for u in "{0}.set_cell = ".format(name)) )

#   # finds with and precision of atoms.
#   if len(getattr(self, atoms)) == 0: return result
#   width = max(len("{0:.0f}".format(u)) for atom in getattr(self, atoms) for u in atom.pos) 
#   precision = ["{0}".format(u) for atom in getattr(self, atoms) for u in atom.pos]
#   precision = max(len(u[(u.index('.') + 1 if '.' in u else 0):]) for u in precision)
#   type_width = max(len(repr(atom.type)) for atom in getattr(self, atoms))
#   start = "{name}.{add_atoms} =".format(name=name, add_atoms=add_atoms)
#   lstart = len(start)
#   site_width, freeze_width = 0, 0
#   for atom in getattr(self, atoms):
#     if atom.site > 0 and len(str(atom.site)) > site_width:
#       site_width = len(str(atom.site))
#     if atom.freeze != FreezeAtom.none and len(str(atom.freeze)) > freeze_width:
#       freeze_width = len(str(atom.freeze))

#   # now prints atoms.
#   for atom in getattr(self, atoms):
#     result += "{start:{lstart}} [({pos[0]:{w}.{p}f}, {pos[1]:{w}.{p}f}, {pos[2]:{w}.{p}f}), "\
#               "{type: >{wtype}}"\
#               .format( start  = start, 
#                        lstart = lstart,
#                        pos    = atom.pos, 
#                        w      = width,
#                        p      = precision,
#                        type   = repr(atom.type),
#                        wtype  = type_width)
#     if atom.site >= 0 or atom.freeze != FreezeAtom.none:
#       if atom.site < 0: result += "None, "
#       else: result += ", {0: >{1}}".format(atom.site, site_width)
#       if atom.freeze != FreezeAtom.none:
#         result += ", {0: >{1}}".format(atom.freeze, freeze_width)
#     result += "],\\\n"
#     start = ""
#   if len(getattr(self, atoms)) > 0: result = result[:-3] + "\n"

#   result += "# End of {name} definition.\n".format(name=name)
#   return result

# def _repr_structure(self): return _repr_lattstr(self)
# _repr_structure.__doc__ = Structure.__str__.__doc__
# Structure.__repr__ = _repr_structure

# def _repr_lattice(self): return _repr_lattstr(self, "sites", "lattice", "add_sites")
# Lattice.__repr__ = _repr_lattice

# def _repr_atom(self): 
#   """ Representation of atoms. """
#   if self.site < 0: 
#     return "{0}({1}, {2})"\
#            .format(self.__class__.__name__, repr(self.pos), repr(self.type))
#   return "{0}({1}, {2}, {3})"\
#          .format(self.__class__.__name__, repr(self.pos), repr(self.type), self.site)
# Atom.__repr__ = _repr_atom
# kAtom.__repr__ = _repr_atom
# rAtom.__repr__ = _repr_atom
# Site.__repr__ = _repr_atom

# def _get_atoms(self): return self._atoms
# def _set_atoms(self, value): self._atoms = self._atoms.__class__(value)
# def _del_atoms(self): self._atoms = self._atoms.__class__([])
# def _get_sites(self): return self._sites
# def _set_sites(self, value): self._sites = self._sites.__class__(value)
# def _del_sites(self): self._atoms = self._atoms.__class__([])
# Lattice.sites = property( _get_sites, _set_sites, _del_sites, \
#                           Lattice._sites.__doc__ )
# Structure.atoms = property( _get_atoms, _set_atoms, _del_atoms, \
#                             Structure._atoms.__doc__ )
# rStructure.atoms = property( _get_atoms, _set_atoms, _del_atoms, \
#                              rStructure._atoms.__doc__ )


# def _copy(self): 
#   """ Returns an exact clone. """
#   from copy import deepcopy
#   return deepcopy(self)

# Lattice.copy   = _copy
# Structure.copy = _copy

# def _copy(self): 
#   """ Returns an exact clone. """
#   from copy import deepcopy
#   return deepcopy(self)
# Lattice.deepcopy   = _copy
# Structure.deepcopy = _copy

# def structure_to_lattice(structure):
#   """ Converts a structure object to a lattice object. """
#   result = Lattice()
#   result.cell = structure.cell.copy()
#   result.name = structure.name
#   result.scale = structure.scale
#   for atom in  structure.atoms:
#     result.add_site = atom.pos, atom.type
#     result.sites[-1].freeze = atom.freeze

#   result.make_primitive()
#   result.find_space_group()
#   # adds  additional components to lattice.
#   result.__dict__.update(structure.__dict__) 
#   return result

# Structure.to_lattice = structure_to_lattice

# def lattice_to_structure(lattice, cell=None, subs=None):
#   """ Creates a structure from a lattice.

#       :Parameters:
#         lattice : `Lattice`
#           Input lattice from which to create structure.
#         cell : None or 3x3 sequence
#           If None, will create structure with primitive cell. Otherwise, will
#           create supercell. In the latter case, the cell should be in cartesian
#           coordinates (not in lattice vector coordinates).
#         subs : None or dict
#           If a dictionary, then substitutes atomic species in the lattice  with
#           other atomic types. E.g. ``subs={'A':'Si'}`` will substitute 'A' in
#           the lattice with 'Si' in the structure. If the lattice site can
#           accomodate more than one atom than the last match will count.
#   """
#   from copy import deepcopy
#   from numpy import array
#   from . import fill_structure
#   if cell  == None: cell = lattice.cell
#   else: cell = array(cell, dtype='float64').reshape((3,3))
#   if subs == None: subs = {}

#   result = fill_structure(cell, lattice)
#   for key, value in subs.items():
#     for atom in result.atoms:
#       if key in lattice.sites[atom.site].type: atom.type = value
#   result.scale = lattice.scale
#   result.name  = lattice.name
#   result.__dict__.update(deepcopy(lattice.__dict__))
#   return result
# Lattice.to_structure = lattice_to_structure


# def vasp_ordered(structure, attributes=None):
#   """ Returns  a structure with correct VASP order of ions.
#   
#       :Parameters: 
#         structure 
#           Structure to reorder
#         attributes
#           A list of attributes to keep in sync with the ions.
#           Defaults to ``["magmom"]``.
#   """
#   from copy import deepcopy
#   if attributes == None: attributes = ["magmom"]
#   elif "magmom" not in attributes and hasattr(attributes, "append"):
#     attributes.append("magmom")

#   result = deepcopy(structure)
#   result.atoms.clear()
#   for attr in attributes:
#     if hasattr(result, attr): setattr(result, attr, [])

#   for type in specie_list(structure):
#     for i, atom in enumerate(structure.atoms):
#       if atom.type != type: continue
#       result.atoms.append(atom)
#       # reorders requested attributes.
#       for attr in attributes:
#         if hasattr(result, attr):
#           getattr(result, attr).append(getattr(structure, attr)[i])
#   return result



# def fill_structure(cell, lattice = None):
#   """ Returns a structure from knowledge of cell and lattice.

#       :Parameters:
#         cell : `Structure`, `rStructure`, or numpy 3x3 float64 array
#           Structure or cell to use to create a complete structure with all atoms.
#         lattice : `Lattice`
#           Back-bone lattice of the super-structure to build. If
#           None, will use the *global* lattice set by `Lattice.set_as_crystal_lattice`.
#       
#       :raises RuntimeError: If the filled structure could not be created.
#   """
#   from _crystal import _fill_structure_impl
#   old_lattice = None
#   try:
#     if lattice == None:
#       try: Structure().lattice
#       except RuntimeError:
#         raise RuntimeError("No lattice given on input of fill_structure" +
#                            "and global lattice not set either.")
#     else: 
#       # keeps track of old lattice.
#       try: old_lattice = Structure().lattice
#       except RuntimeError: pass
#       # sets this lattice as the global lattice.
#       lattice.set_as_crystal_lattice()

#     # creates filled structure.
#     result = _fill_structure_impl(cell)
#   except: raise
#   finally: 
#     # Now resets old lattice.
#     if old_lattice != None: old_lattice.set_as_crystal_lattice()
#   
#   return result

# Structure.write_poscar = write_poscar
# Structure.write_oldvff = write_oldvff

# def gaussian_projector(positions, cell, center=(0,0,0), alpha=1e0):
#   """ Returns gaussian projector operator. 
#   
#       :Parameters:
#         positions : numpy array
#           Last dimension of this array must be 3. 
#         center : numpy array
#           1-dimensional vector of three elements. Center of the gaussian projector.
#         cell : numpy array
#           Cell into which gaussian projector is wrapped. I.e. this operator is
#           periodic.
#         alpha : numpy scalar
#           Factor within the gaussian, akin (but not equal) to width at half maximum.

#       :return: array = exp( alpha * norm(positions[i])**2 )

#       All arguments can have units. The return however is without units
#       (dimensionless).
#   """
#   from _crystal import _gaussian_projector_impl
#   from numpy import array
#   length_units = None
#   if isinstance(center, tuple) or isinstance(center, list): center = array(center)
#   if hasattr(positions, "units"): length_units = positions.units
#   elif hasattr(center, "units"):  length_units = center.units
#   elif hasattr(cell, "units"):    length_units = cell.units
#   if length_units != None:
#     if hasattr(cell, "units"):
#       try: cell = cell.rescale(length_units)
#       except: raise ValueError("Cell and positions have incompatible units.")
#     if hasattr(center, "units"):
#       try: center = center.rescale(length_units)
#       except: raise ValueError("Center and positions have incompatible units.")
#     if hasattr(alpha, "units"):
#       try: alpha_units = length_units ** -2
#       except: raise ValueError("Could not create alpha units from positions units.")
#       try: alpha = alpha.rescale(alpha_units)
#       except: raise ValueError("Could not rescale alpha to {0}.".format(alpha_units))
#       else: alpha = float(alpha.magnitude)
#   return _gaussian_projector_impl(positions, cell, center, alpha)

# @contextmanager
# def lattice_context(lattice):
#   """ Sets global lattice within a context. """
#   from _crystal import _nullify_global_lattice
#   from . import Structure
#   try: oldlattice = Structure().lattice
#   except: oldlattice = None
#   lattice.set_as_crystal_lattice()
#   yield oldlattice
#   if oldlattice != None: oldlattice.set_as_crystal_lattice()
#   else: _nullify_global_lattice()

# def layer_iterator(structure, direction, tolerance=1e-12):
#   """ Iterates over layers and atoms in a layer. 

#       :Parameters: 
#         structure : `Structure`
#           The structure for which to iterator over atoms.
#         direction : 3d-vector
#           Growth direction vector, e.g. vector perpendicular to the layers.
#           Defaults to the first column vector of the structure.  It is
#           important that two cell-vectors of the structure are (or can be
#           transformed to be) perpendicular to the growth direction. Otherwise
#           it cannot be assured that layers are well defined, i.e. that each
#           atom belongs to a single (as in periodic) layer. This condition is
#           *not* enforced (no assertion) since it is not necessary, only
#           sufficient.  Note that the third vector need not be parallel to the
#           growth direction.
#         tolerance : float
#           Maximum difference between two atoms in the same layer.
#   """
#   from operator import itemgetter
#   from numpy import array, dot
#   from . import LayerDepth, to_cell, to_voronoi

#   direction = array(direction)
#   if len(structure.atoms) <= 1: yield structure.atoms; return

#   # orders position with respect to direction.
#   positions = to_cell(array([atom.pos for atom in structure.atoms]), structure.cell)
#   projs = [(i, dot(pos, direction)) for i, pos in enumerate(positions)]
#   projs = sorted(projs, key=itemgetter(1))

#   # creates classes of positions.
#   result = [[projs[0]]]
#   for i, proj in projs[1:]:
#     if abs(proj - result[-1][-1][-1]) < tolerance: result[-1].append((i,proj))
#     else: result.append([(i,proj)])

#   # only one layer.
#   if len(result) == 1: yield structure.atoms; return
#   # Finds if first and last have atoms in common through periodicity
#   first, last = result[0], result[-1]
#   centered = to_voronoi(positions[[i for i, d in last]], structure.cell, positions[first[0][0]])
#   for j, pos in enumerate(centered[::-1]):
#     a0 = dot(pos, direction)
#     if any(abs(u[1]-a0) >= tolerance for u in first): continue
#     first.append( last.pop(len(centered)-j-1) )

#   # last layer got depleted.
#   if len(last) == 0: result.pop(-1) 
#   # only one layer.
#   if len(result) == 1: yield structure.atoms; return
#   # yield layer iterators.
#   for layer in result:
#     def inner_layer_iterator():
#       """ Iterates over atoms in a single layer. """
#       for index, norm in layer: yield structure.atoms[index]
#     yield inner_layer_iterator ()


# def equivalence_iterator(structure, operations = None, tolerance=1e-6):
#   """ Yields iterators over equivalent atoms.
#   
#       :Parameters:
#         structure : `Structure` or `Lattice`
#           Structure or Lattice over which to iterate.
#         operations : Iterable over symmetry operations
#           A symmetry operation is a callable which takes a 3d-vector and
#           returns a 3d-vector. Equivalence is judged according to positions
#           only. In general, this parameter will be ``lattice.space_group``.
#           If None, the operations will be set as follows:
#           
#             - structure is a `Structure` instance: a lattice is created from
#               the structure, and its space group operations are used.
#             - structure is a `Lattice` instance: its space group operations are used.

#         tolerance : float
#           Two positions closer than ``tolerance`` are considered equivalent.
#   """
#   from numpy import array
#   from numpy.linalg import norm
#   from . import to_cell

#   if hasattr(structure, 'to_lattice'):
#     atoms = [u for u in enumerate(structure.atoms)]
#     is_structure = True
#     if operations == None: operations = structure.to_lattice().space_group
#   else:
#     atoms = [u for u in enumerate(structure.sites)]
#     is_structure = False
#     if operations == None: operations = structure.space_group
#    
#   while len(atoms):
#     i, atom = atoms.pop()
#     equivs = [i]
#     if len(atoms): 
#       for op in operations:
#         others = to_cell(array([u[1].pos for u in atoms]), structure.cell, op(atom.pos))
#         others = [i for i, pos in enumerate(others) if norm(pos) < tolerance]
#         for index in others:
#           i, pos = atoms.pop(index)
#           equivs.append(i)
#     if is_structure:
#       def inner_equivalence_iterator():
#         """ Iterates over equivalent atoms in a structure. """
#         for i in equivs: yield structure.atoms[i]
#     else:
#       def inner_equivalence_iterator():
#         """ Iterates over equivalent atoms in a structure. """
#         for i in equivs: yield structure.sites[i]
#     yield inner_equivalence_iterator()
