""" Contains basic data type and methods for crystal structures. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Structure', 'Atom', 'HFTransform', 'zero_centered', 'into_voronoi',
            'into_cell', 'supercell', 'primitive', 'is_primitive', 'space_group',
            'transform', 'periodic_dnc', 'neighbors', 'coordination_shells',
            'splitconfigs', 'vasp_ordered', 'specieset', 'map_sites',
            'which_site', 'iterator' ]
from cppwrappers import Structure, Atom, HFTransform, zero_centered, into_voronoi,    \
                        into_cell, supercell, primitive, is_primitive, space_group,   \
                        transform, periodic_dnc, neighbors, coordination_shells,      \
                        splitconfigs, map_sites
import iterator

def specieset(structure):
  """ Returns ordered set of species.
  
      Especially usefull with VASP since we are sure what the list of species
      is always ordered the same way.
  """
  return set([a.type for a in structure])

def vasp_ordered(structure):
  """ Returns a structure with correct VASP order of ions.
  
      :param structure:
          :class:`Structure` for which to reorder atoms.
  """

  from copy import deepcopy
  result = deepcopy(structure)
  def sortme(self): return a.type.lower()
  result[:] = sorted(structure, sortme)
  return result

def which_site(atom, lattice, invcell=None, tolerance=1e-8):
  """ Index of periodically equivalent atom. 


      :param atom: 
        :py:class:`~cppwrappers.Atom` for which to find periodic equivalent.
      :param lattice:
        :py:class:`~cppwrappers.Structure` defining the periodicity.
      :type lattice: :py:class:`~cppwrappers.Structure` or matrix

      :return: index in list of atoms, or -1 if not found.
  """
  from numpy.linalg import inv
  from .cppwrappers import are_periodic_images as api
  if invcell is None: invcell = inv(lattice.cell)
  lattice = [getattr(site, 'pos', site) for site in lattice]
  pos = getattr(atom, 'pos', atom)
  for i, site in enumerate(lattice):
    if api(pos, site, invcell, tolerance): return i
  return -1
