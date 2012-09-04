""" Contains basic data type and methods for crystal structures. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Structure', 'Atom', 'HFTransform', 'zero_centered', 'into_voronoi',
            'into_cell', 'supercell', 'primitive', 'is_primitive', 'space_group',
            'transform', 'periodic_dnc', 'neighbors', 
            'coordination_shells', 'splitconfigs', 'vasp_ordered', 'layer_iterator',
            'equivalence_iterator', 'shell_iterator', 'specieset', 'map_sites',
            'which_site' ]
from cppwrappers import Structure, Atom, HFTransform, zero_centered, into_voronoi,    \
                        into_cell, supercell, primitive, is_primitive, space_group,   \
                        transform, periodic_dnc, neighbors, coordination_shells,      \
                        splitconfigs, map_sites

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

def layer_iterator(structure, direction, tolerance=1e-12):
  """ Iterates over layers and atoms in a layer. 

      :param structure: 
          :class:`Structure` for which to iterator over atoms.
      :param direction:
          3d-vector defining the growth direction, e.g. vector perpendicular to the layers.
          Defaults to the first column vector of the structure.  It is
          important that two cell-vectors of the structure are (or can be
          transformed to be) perpendicular to the growth direction. Otherwise
          it cannot be assured that layers are well defined, i.e. that each
          atom belongs to a single (as in periodic) layer. This condition is
          *not* enforced (no assertion) since it is not necessary, only
          sufficient.  Note that the third vector need not be parallel to the
          growth direction.
      :param tolerance: 
          Maximum difference between two atoms in the same layer.

      :returns: Yields iterators over atoms in a single epitaxial layer.
  """
  from operator import itemgetter
  from numpy import array, dot

  direction = array(direction)
  if len(structure) <= 1: yield list(structure); return

  # orders position with respect to direction.
  positions = array([into_cell(atom.pos, structure.cell) for atom in structure])
  projs = [(i, dot(pos, direction)) for i, pos in enumerate(positions)]
  projs = sorted(projs, key=itemgetter(1))

  # creates classes of positions.
  result = [[projs[0]]]
  for i, proj in projs[1:]:
    if abs(proj - result[-1][-1][-1]) < tolerance: result[-1].append((i,proj))
    else: result.append([(i,proj)])

  # only one layer.
  if len(result) == 1: yield structure; return
  # Finds if first and last have atoms in common through periodicity
  first, last = result[0], result[-1]
  centered                                                                     \
    = into_voronoi( positions[[i for i, d in last]] - positions[first[0][0]],
                    structure.cell )
  for j, pos in enumerate(centered[::-1]):
    a0 = dot(pos, direction)
    if any(abs(u[1]-a0) >= tolerance for u in first): continue
    first.append( last.pop(len(centered)-j-1) )

  # last layer got depleted.
  if len(last) == 0: result.pop(-1) 
  # only one layer.
  if len(result) == 1: yield structure; return
  # yield layer iterators.
  for layer in result:
    def inner_layer_iterator():
      """ Iterates over atoms in a single layer. """
      for index, norm in layer: yield structure[index]
    yield inner_layer_iterator ()


def equivalence_iterator( structure, operations=None,                          \
                          tolerance=1e-6, splitocc=None ):
  """ Yields iterators over atoms equivalent via space group operations.
  
      Only check that the position are equivalent. Does not check that the
      occupations are the same.

      :param structure:
          :class:`Structure` over which to iterate.
      :param operations: 
          List of symmetry operations.
          A symmetry operation is 4x4 matrix where the upper block is a
          rotation and the lower block a translation. The translation is
          applied *after* the rotation. If None, the operations are obtained
          using:class:`space_group`.
      :param float tolerance:
          Two positions closer than ``tolerance`` are considered equivalent.
      :param callable splitocc:
          Function to split two sites according to something other than
          geometry. Generally, this would be occupation and/or magnetic state.
          It is a callable taking two :py:class:`cppwrappers.Atom` and
          returning True or False depending on whether they are equivalent.
          It should be transitive and symmetric, otherwise results will be
          undetermined.  If None, then splitting occurs only according to
          geometry.
      
      :returns: Yields iterators over atoms linked by space-group operations.
  """
  from numpy import array, dot
  from numpy.linalg import norm
  from . import into_cell

  atoms = [u for u in enumerate(structure)]
  if operations == None: operations = space_group(structure)
   
  while len(atoms):
    i, atom = atoms.pop()
    equivs = [i]
    if len(atoms): 
      for op in operations:
        others = into_cell(array([u[1].pos for u in atoms]),                   \
                 structure.cell, dot(op[:3], atom.pos)) + op[3]
        others = [i for i, pos in enumerate(others) if norm(pos) < tolerance]
        for index in others:
          i, pos = atoms.pop(index)
          equivs.append(i)
    if splitocc == None: yield equivs
    else: 
      results = []
      for u in equivs: 
        if len(results) == 0: results.append(equivs)
        else:
          found = False
          for group in results:
            if splitocc(group[0], u):
              group.append(u)
              found = True
              break
          if not found: results.append([u])
      for u in results: yield u

def shell_iterator(structure, center, direction, thickness=0.05):
  """ Iterates over cylindrical shells of atoms.
  
      It allows to rapidly create core-shell nanowires.
  
      :param structure: 
          :class:`Structure` over which to iterate.
      :param center: 
          3d-vector defining the growth direction of the nanowire.
      :param thickness: 
          Thickness in units of ``structure.scale`` of an individual shell.
      
      :returns: Yields iterators over atoms in a single shell.
  """
  from operator import itemgetter
  from numpy import array, dot
  from numpy.linalg import norm
  from operator import into_voronoi

  direction = array(direction)/norm(array(direction))
  if len(structure) <= 1: yield structure; return

  # orders position with respect to cylindrical coordinates
  positions = into_voronoi(array([atom.pos - center for atom in structure]), structure.cell)
  projs = [(i, norm(pos - dot(pos, direction)*direction)) for i, pos in enumerate(positions)]
  projs = sorted(projs, key=itemgetter(1))

  # creates classes of positions.
  result = {}
  for i, r in projs:
    index = int(r/thickness+1e-12)
    if index in result: result[index].append(i)
    else: result[index] = [i]

  for key, layer in sorted(result.iteritems(), key=itemgetter(0)):
    def inner_layer_iterator():
      """ Iterates over atoms in a single layer. """
      for index in layer: yield structure[index]
    yield inner_layer_iterator()

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
  pos = getattr(atom, 'pos', atom)
  for i, site in enumerate(lattice):
    if api(pos, site.pos, invcell, tolerance): return i
  return -1
