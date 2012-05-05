def objective(self, extractor, indiv):
  """ Saves info into invididual to later compute distance from convex-hull. """
  from numpy import array
  from quantities import eV
  stoic = [len([0 for a in indiv if a.type == t]) for t in self.species]
  indiv.concentrations = array(stoic) / float(len(indiv))
  if hasattr(extractor.total_energy, 'rescale'): 
    indiv.energy = float(extractor.total_energy.rescale(eV).magnitude)
  else: indiv.energy = extractor.total_energy
  return indiv.energy


def strip_points(points):
  """ Strip away points at same concentration, keeping only lowest energy. """ 
  from numpy import array, all, abs
  result = [points[0]]
  for point in points[1:]:
    found = False
    for i, cmp in enumerate(result):
      if all(abs(array(point[:-1]) - cmp[:-1]) < 1e-8):
        found = True
        if point[-1] < cmp[:-1]: result[i] = point
        break
    if not found: result.append(point)
  return array(result)


def update_fitness(self, individuals):
  """ Returns distance from convex-hull.

      :param ch:
        Convex-hull as a Vrep from Polyhedron_ python package.
      :param x:
        Concentration coordinates
      :param y:
        Energy coordinate.
      
      .. _Polyhedron:: http://cens.ioc.ee/projects/polyhedron
  """
  from itertools import chain
  from numpy import array, concatenate, dot, max
  from polyhedron import Vrep
  individuals = list(individuals)
  # coordinates of individuals
  points = [u.concentrations[:-1].tolist() + [u.energy] for u in individuals]
  # if a convex-hull already exists, adds generators.
  if hasattr(self, 'convexhull'): 
    allpoints = set(chain(self.convexhull.ininc))
    points.extend(self.generators[list(allpoints)])

  if len(points) == 0: return

  # strip away points at same stiochiometry, keeping only lowest.
  points = strip_points(points)

  # create convex hull from all points.
  vrep = Vrep(points)
  
  # lowers are those half-spaces on the bottom of the energy scale.
  lowers = concatenate((vrep.A, vrep.b[:,None]), axis=1)[vrep.A[:,-1].flat<-0.1]

  # now update individuals
  for indiv in individuals:
    point = array(indiv.concentrations[:-1].tolist() + [indiv.energy])
    indiv.fitness = max(point[-1] -(dot(lowers[:,:-2], point[:-1]) - lowers[:,-1] ) / lowers[:,-2])

def random_cell(angle_range=None, cubic=False, **kwargs):
  """ Generates random cell. 
      
      In the algorithm below, a,b,c are the lengths of the lattice vectors
      alpha, beta, gamma are the angles betweeen them - the standard convetnion
      is respected (see `here <http://en.wikipedia.org/wiki/Crystal_structure>`_)
      lengths are between 0.6 and 1.4 in units of scale angles in the
      angle_range range

      :param float angle_range: 
           Minimum and maximum angles between lattice vectors.
           Defaults to 60.0 and 140.0 degrees.
      :param bool cubic: 
           Whether the cell should be cubic or not.
      :param kwargs:
           Passed on to :py:func:`~lada.math.gruber`

      :returns: A random crystal cell.
  """
  from random import random
  from numpy import pi, sqrt, cos, sin, array, zeros, abs
  from numpy.linalg import det
  from ..math import gruber

  if angle_range == None: angle_range = [60.,140.]
  cell = zeros([3,3])
  while abs(det(cell)) < 1e-6:

    if cubic:
      a = b = c = 1.
      alpha = beta = gamma = pi/2.
    else:
      a = 0.8*random()+0.6 
      b = 0.8*random()+0.6
      c = 0.8*random()+0.6
        
      alpha = ((angle_range[1]-angle_range[0])*random() + angle_range[0])*pi/180.
      beta  = ((angle_range[1]-angle_range[0])*random() + angle_range[0])*pi/180.
      gamma = ((angle_range[1]-angle_range[0])*random() + angle_range[0])*pi/180.

    a1 = a*array([1.,0.,0.])
    a2 = b*array([cos(gamma),sin(gamma),0.])
        
    c1 = c*cos(beta)
    c2 = c/sin(gamma)*(-cos(beta)*cos(gamma) + cos(alpha))
        
    a3 = array([c1, c2, sqrt(c**2-(c1**2+c2**2))])
        
    cell = array([a1,a2,a3]).T
    if det(cell) < 0e0: cell = array([a2, a1, a3]).T
  return gruber(cell.copy(), **kwargs)

def random_structure( nbatoms, min_distance=1.3, atomic_scale=None,
                      specie='anion', **kwargs): 
  """ Generates random crystal structure.
      
      Creates a random cell using :py:meth:`random_cell`, then places atoms
      randomly inside the cell.
      contraints on the range of angles between lattice 
      vectors and the minimal distance distances between 
      two atoms.
      
      :param int nbatoms: 
         Number of atoms in the structure.
      :param real min_distance: 
         Minimal distance between two atoms. If no units, defaults to angstrom.
      :param real atomic_scale:
         Used in probability grid when an atom is placed. A guassian with this
         half-width is created which makes placing an atom here less probable.
         If None, then we use min_distance.
      :param **kwargs: 
         Extra parameters are passed on to :py:meth:`random_cell`.

      :returns: A random crystal structure.
  """
  from random import choice
  from numpy import dot, array, arange, ones, argwhere, exp, sum, argmin
  from itertools import product, chain
  from quantities import angstrom
  from lada.crystal import Structure, neighbors

  # creates structure with random cell.
  result = Structure(random_cell(**kwargs))
  # sets scale using min_distance as a first guess.
  if hasattr(min_distance, 'rescale'):
    min_distance = float(min_distance.rescale(angstrom))
  minvol = 4.0/3.0*min_distance**3 * nbatoms
  result.scale = (2*minvol/result.volume)**(1.0/3.0)
  # set atomic_scale to result.scale units.
  if hasattr(atomic_scale, 'rescale'): 
    atomic_scale = float(atomic_scale.rescale(angstrom))
  if atomic_scale is None: atomic_scale = min_distance
  atomic_scale /= result.scale

  # Creates a grid where atoms will be placed.
  N = max(7, (3.0*nbatoms**(1./3.)))
  grid = array(list(product( arange(0, 1., 1.0/float(N)), 
                             arange(0, 1., 1.0/float(N)), 
                             arange(0, 1., 1.0/float(N)) )))
  grid = dot(result.cell, grid.T).T
  
  # Setting the initial prob. distribution
  probability  = ones(len(grid), dtype="float64")
  
  # periodic images. used to construct periodict probability density where
  # atoms shouldn't be.
  periodic_images = [ dot(result.cell, [i,j,k]) for i in arange(-1,2)\
                      for j in arange(-1,2) for k in arange(-1,2) ]
  periodic_images = array(periodic_images)

  # places first atom at origin.
  result.add_atom(0, 0, 0, type=specie)
  # make neighbors less probable
  distances = grid[:,None, :] - periodic_images[None, :, :] 
  probability -= sum(exp(-sum(distances**2, axis=2)/atomic_scale**2), axis=1)

  # now adds each atom in turn. 
  for i in xrange(1, nbatoms):
    # create criteria for placing atoms.
    maxi, mini = max(probability), min(probability)
    criteria = 0.90 * (maxi - mini) + mini
    # Find points on the anion grid that have the desired probability
    empty_space=argwhere(probability >= criteria) 
    # choose position randomly
    pom = choice(empty_space.flatten())
    # add atom to structure
    result.add_atom(*grid[pom].flatten(), type=specie)
    # make neighbors less probable
    distances = (grid- grid[pom])[:,None, :] - periodic_images[None, :, :] 
    probability -= sum(exp(-sum(distances**2, axis=2)/atomic_scale**2), axis=1)

  # Set the scale so that the min_distance is fulfilled
  mindist = [ (vector, distance) for atom in result \
              for neighbor, vector, distance in neighbors(result, 1, atom) ]
  mindist = mindist[argmin([d[1] for d in mindist])]
  if mindist[1] * result.scale < min_distance:
    result.scale *= min_distance / result.scale / mindist[1] 

  return result

def populate_anion_cation(structure, species, anions):
  """ Populates structure with anions and cations. 
  
      Tries to place cations in such a way as to break as many anion-anion
      bonds. Then places the anion and cation species randomly within each
      sublattice.
  
      :param structure:
        :py:class:`~lada.crystal.Structure` to modify. All atomic types are
        replaced.
      :param dict species: 
        Dictionary where keys correspond to atomic species, and values to the
        stoichiometry.
      :param dict anions:
        List of anionic species. Cations are the keys in ``species`` not listed
        in ``anions``.
  """
  from random import shuffle, choice
  from numpy import array, sum, argmax
  from ..crystal import coordination_shells

  # first creates list of cations and anions.
  types = anions
  anions, cations = [], []
  for key, value in species.iteritems():
    if key not in types: 
      cations.extend([key]*value)
    else: anions.extend([key]*value)
  if len(anions) + len(cations) != len(structure): 
    raise ValueError( "Number of species and atoms do not match.\n"\
                      "There are {0} anions and {1} cations, but "\
                      "{2} atoms."\
                      .format(len(anions), len(cations), len(structure)) )
    
  # makes sure all atoms are called anions.
  if len(anions) == 0: # unless no anions are requested.
    for a in structure: a.type = 'cation'
  else: 
    for atom in structure: atom.type = 'anions'
  # if cations exist, place them so as to minimize the number of cation-cation
  # and anion-anion bonds.
  if len(cations) != 0:
    # loop over shells, to minimize same-same bonds.
    shells = [ coordination_shells(structure, 1, atom, 0.5/structure.scale)[0] \
               for atom in structure ]
    # Each iteration places another cation.
    for i in xrange(len(cations)):
      # environment holds the same-same bond score.
      # since that score changes each time a cation is placed, it is
      # recomputed.
      environment = []
      for shell, atom in zip(shells, structure):
        if atom.type == 'cation': x = -1 # avoids double placement.
        else: 
          # score if function of bond-length.
          x = [(1.0/d if atom.type == n.type else 0.) for n, v, d in shell]
        environment.append(sum(x))
      environment = array(environment)
      # finds all highest scoring environments, and chooses a highest-scorer
      # randomly.
      indices = environment.argsort(axis=0)
      u = argmax(environment[indices])
      u = choice(indices[u:])
      structure[u].type = 'cation'

  # now replaces 'anion'/'cation' with actual species.
  shuffle(anions)
  shuffle(cations)
  for atom in structure:
    if atom.type == 'cation':
          atom.type = cations.pop()
    else: atom.type = anions.pop()

def populate(structure, species):
  """ Randomly populates structure with species. """
  from random import shuffle
  types = []
  for key, value in species.iteritems(): types.extend([key]*value)
  for atom, type in zip(structure, types): atom.type = type

def taboo( structure, max_atoms=-1, min_distance=1.3, \
           angle_range=None, same_first_neigh=2, verbose=True ):
  """ Checkes geometry of the structure and returnes 

      True if the structure satisfyies a set of conditions
      otherwise returns False
      
      :param structure:
          :py:class:`~lada.crystal.Structure` to check.
      :param int max_atoms: 
          Maximal allowed number of atoms
      :param float min_distance: 
          Minimal first shell distance in angstroms
      :param angle_range: 
          Range of allowed angles between lattice vectors
      :type angle_range: 2x1 float list
      :param int same_first_neigh: 
          Maximal number of atoms that have more than two 
          neighbors of the same type whithin 1.5*min_distance.
          If strictly negative, then does not check.
      :param bool verbose: 
          Whether to print messages
  """
  from numpy import pi, arccos, dot, sqrt
  from ..crystal import Neighbors
  
  # check number of atoms.
  if max_atoms > 0 and len(structure) > max_atoms: 
    if verbose: print "more than ",max_atoms," atoms"
    return True

  # chekc first neighbor distances.
  distances = []
  for atom in structure:
    for n in Neighbors(structure,1,atom.pos):
      distances.append(n.distance)

  if structure.scale*min(distances) < min_distance: 
    if verbose: print "distances shorter than",min_distance," AA"
    return True

  # check that angles are within specified range.
  if angle_range is None: angle_range = (45.,160.)
  cell = structure.cell
  alpha = 180./pi*arccos( dot(cell[:,0],cell[:,2])/sqrt(dot(cell[:,0],cell[:,0])*dot(cell[:,2],cell[:,2])) )
  beta  = 180./pi*arccos( dot(cell[:,1],cell[:,2])/sqrt(dot(cell[:,1],cell[:,1])*dot(cell[:,2],cell[:,2])) )
  gamma = 180./pi*arccos( dot(cell[:,0],cell[:,1])/sqrt(dot(cell[:,0],cell[:,0])*dot(cell[:,1],cell[:,1])) )
  
  if not angle_range[0] <= alpha <= angle_range[1]:
    if verbose: print "alpha =",alpha," degrees"
    return True

  if not angle_range[0] <= beta  <= angle_range[1]:
    if verbose: print "beta =",beta," degrees"
    return True

  if not angle_range[0] <= gamma <= angle_range[1]:
    if verbose: print "gamma =",gamma," degrees"
    return True


  # check neighbor types.
  if same_first_neigh >= 0:
    check = []
    for atom in structure:
      pom = 0
      for n in Neighbors(structure,20,atom.pos):
        if structure.scale*n.distance <= 1.5*min_distance and structure[n.index].type==atom.type:
          pom = pom + 1
      
      check.append(pom)

    how_many = len([x for x in check if x>2])

    if how_many > same_first_neigh: 
      if verbose: print how_many," atoms with more than 2 close neighbors of the same type"
      return True

  # all done
  return False

def fractional_pos(structure, roll=True):
  """ Returns fractional coordinates of atoms in structure. """
  from random import random
  from numpy import dot, array
  from numpy.linalg import inv
  # random translation vector.
  trans = (random(), random(), random()) if roll else (0,0,0)
  result = [dot(inv(structure.cell), atom.pos) + trans for atom in structure]

  for pos in result: pos = mod(pos, [1,1,1])

  return array(result)

def cut_and_splice(s1, s2, roll=True):
  """ Cut-n-splice GSGO crossover operation

      Mating operation on two parent structures s1 and s2.
      Done on scaled atomic positions (cubic systems) by
      cutting them in half and mixing their upper and 
      lower parts.
  """
  from random import choice, random
  from numpy import dot, abs
  from numpy.linalg import det
  from lada.crystal import Structure

  # swap structures randomly.
  if choice([True, False]): s1, s2 = s2, s1

  # chem. symbols and scaled positions of the two parents
  sc_pos1 = zip([atom.type for atom in s1],fractional_pos(s1, roll=True))
  sc_pos2 = zip([atom.type for atom in s2],fractional_pos(s2, roll=True))

  result = Structure(s1.cell, scale=s1.scale)

  # choose random positions of split-plane
  xsep = 0.5 - (random() * 0.45 + 0.15)
  # choose direction of split-plane randomly from cell-vectors.
  direction = choice(range(3))

  for type, pos in sc_pos1:
    if pos[direction] >= xsep: result.add_atom(*dot(result.cell, pos), type=type)

  for type, pos in sc_pos2:
    if pos[direction] < xsep: result.add_atom(*dot(result.cell, pos), type=type)

  result.scale =  s1.scale * (float(len(result)) / float(len(s1)))**(1./3.)

  return result


def mix_poscars(s1,s2, roll=True):
  """ Crossover operations where atoms are from one and cell from other.

      Mating operation on two parent structures s1 and s2.
      Done on scaled atomic positions (cubic systems) by
      interchanging their cells and scaled positions.
      Returns two offspring structures each with the same 
      number of atoms as the parent from which the atoms 
      are inhereted.
  """
  from random import choice
  from numpy import dot
  from numpy.linalg import det
  from ..crystal import Structure

  # swap structures randomly.
  if choice([True, False]): s1, s2 = s2, s1

  # chem. symbols and scaled positions of the two parents
  sc_pos2 = zip([atom.type for atom in s2],fractional_pos(s2, roll))

  # cell from s1
  result = Structure(s1.cell, scal=s1.scale)

  # atoms from s2
  for type, pos in sc_pos2:
    result.add_atom(*dot(result.cell, pos), type=type)

  result.scale =  s1.scale * (float(len(result)) / float(len(s1)))**(1./3.)
  return result

def mix_atoms(s1, s2, roll=True):
  """ Randomly mix cell and atoms from parents.

      Mating operation on two parent structures s1 and s2.
      Done by mixing randomly atoms from s1 and s2 into the 
      cell coming from one of them. Returns two offspring 
      structures.
  """
  from random import choice
  from itertools import chain
  from numpy.linalg import det
  from numpy import dot, abs
  from ..crystal import Structure

  # swap structures randomly.
  if choice([True, False]): s1, s2 = s2, s1

  # chem. symbols and scaled positions of the two parents
  sc_pos1 = zip([atom.type for atom in s1], fractional_pos(s1, roll))
  sc_pos2 = zip([atom.type for atom in s2], fractional_pos(s2, roll))

  result = Structure(s1.cell)

  for pos, type in chain(sc_pos1, sc_pos2):
    if choice([True, False]):
      result.add_atom(*dot(result.cell, pot), type=type)

  result.scale =  s1.scale * (float(len(result)) / float(len(s1)))**(1./3.)

  return result

def jiggle_structure(structure):
  """ Mutates by jiggling structure. """
  from copy import deepcopy
  from numpy import dot, min
  from numpy.linalg import norm, det
  from numpy.random import rand
  result = deepcopy(structure)
  omega = rand(3,3) * 2e0 - 1e0
  omega += omega.T
  newcell = result.cell + dot(omega, result.cell)

  fractionals = fractional_pos(result)
  mindists = [min([norm(a.pos - atom.pos) for a in result]) for atom in result]
  for atom, frac, mindist  in zip(result, fractionals, mindists):
    atom.pos = dot(newcell, frac)  + mindist * 0.25 * (rand(3) * 2e0 - 1e0)
  result.cell = newcell

  result.scale = structure.scale * (abs(det(structure.cell)) / abs(det(result.cell))\
                  * float(len(result)) / float(len(structure)))**(1./3.)
  return result
