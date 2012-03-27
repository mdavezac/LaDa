def objective(self, extractor, indiv):
  """ Saves info into invididual to later compute distance from convex-hull. """
  from numpy import array
  from quantities import eV
  stoic = [len([0 for a in indiv.atoms if a.type == t]) for t in self.species]
  indiv.concentrations = array(stoic) / float(len(indiv.atoms))
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
  from numpy import array, concatenate, dot
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
    print point[-1] -(dot(lowers[:,:-2], point[:-1]) - lowers[:,-1] ) / lowers[:,-2] 
    indiv.fitness = point[-1] -(dot(lowers[:,:-2], point[:-1]) - lowers[:,-1] ) / lowers[:,-2] 

def rand_cell_gen(angles_between=None, cubic=False):
  """ Returns random cell. 
      
      In the algorithm below, a,b,c are the lengths of the lattice vectors
      alpha, beta, gamma are the angles betweeen them - the standard convetnion
      is respected (see `here <http://en.wikipedia.org/wiki/Crystal_structure>`_)
      lengths are between 0.6 and 1.4 in units of scale angles in the
      angles_between range

      :param angles_between: 
           Minimum and maximum angles between lattice vectors.
           Defaults to 60.0 and 140.0

      :returns: A random crystal structure.
  """
  from random import random
  from numpy import pi, sqrt, cos, sin, array, zeros, abs
  from numpy.linalg import det

  if angles_between == None: angles_between = [60.,140.]
  cell = zeros([3,3])
  while abs(det(cell)) < 1e-6:

    if cubic:
      a = b = c = 1.
      alpha = beta = gamma = pi/2.
    else:
      a = 0.8*random()+0.6 
      b = 0.8*random()+0.6
      c = 0.8*random()+0.6
        
      alpha = ((angles_between[1]-angles_between[0])*random() + angles_between[0])*pi/180.
      beta  = ((angles_between[1]-angles_between[0])*random() + angles_between[0])*pi/180.
      gamma = ((angles_between[1]-angles_between[0])*random() + angles_between[0])*pi/180.

    a1 = a*array([1.,0.,0.])
    a2 = b*array([cos(gamma),sin(gamma),0.])
        
    c1 = c*cos(beta)
    c2 = c/sin(gamma)*(-cos(beta)*cos(gamma) + cos(alpha))
        
    a3 = array([c1, c2, sqrt(c**2-(c1**2+c2**2))])
        
    cell = array([a1,a2,a3])
  if det(cell) < 0e0: cell = array([a2, a1, a3])
  return cell.T

def rand_struct_gen( species, stoichiometry, anions=None, angles_between=None,
                     min_distance=1.7, atomic_scale=2., cubic=False):
  """ Function that returnes random structure within 
      
      contraints on the range of angles between lattice 
      vectors and the minimal distance distances between 
      two atoms.
      
      :param species: 
           Chemical symbols of the atoms involved
      :type species: list of n strings
      :param stoichiometry: nx1 int list 
           Number of atoms of each type
      :type stoichiometry: list of n strings
      :param anions: 
           List of zeros and ones denoting cations and anions, respectively.
           Ignored if None
      :param angles_between: 
           Minimum and maximum angles between lattice vectors.
           Defaults to 60.0 and 140.0
      :param real min_distance: 
           Minimal distance bwtween two atoms
      :param real atomic_scale:
           Factors the minimal grid distance to determine the ionic size 
           Be aware this has nothing to do with real atomic/ionic sizes
           Change it if the atoms appear to be too close

      :returns: A random crystal structure.
  """
  from random import randint,choice
  from lada.crystal import Structure, Neighbors
  from numpy import pi, sqrt, cos, dot, array, transpose, arange
  from numpy.linalg import inv

  def prob(x,width,center,structure):
    from numpy import exp, dot, array, arange
    centers = [ center + dot(structure.cell, [i,j,k]) for i in arange(-1,2)\
                for j in arange(-1,2) for k in arange(-1,2) ]
    return sum(array([exp(-dot(x-c,x-c)/width**2) for c in centers]))

  
  if anions == None: anions = [1]*len(species)
  result = Structure()
  result.scale = 1.

  cell = rand_cell_gen(angles_between, cubic)
  rec_cell = 2*pi*transpose(inv(cell)) # reciprocal cell

  result.cell=transpose(cell)

  # Setting two grids, one for the anions and one for the cations
  # This is done by taking the maxima and minima points of the COS(g*r) function 
  # where g is a random reciproical lattice vector (in other words a supperlattice)
  grid_an =[]
  grid_cat=[]

  # setting the g
  # chioce between 4 and 7 is to have large enough number of points, but not too large grids
  slat = [randint(4,7),randint(4,7),randint(4,7)]
  g = slat[0]*rec_cell[0]+slat[1]*rec_cell[1]+slat[2]*rec_cell[2]

  # setting the real-space grids in such a way that they contain all the mentioned maxima and minima points
  grid_x=[x for x in arange(0.,1.,1./(2*slat[0]))]
  grid_y=[x for x in arange(0.,1.,1./(2*slat[1]))]
  grid_z=[x for x in arange(0.,1.,1./(2*slat[2]))]

  # Filling the candidate points for anions (maxima) and cations (minima)
  for i in range(len(grid_x)):
    for j in range(len(grid_y)):
      for k in range(len(grid_z)):
        pom=dot(result.cell,array([grid_x[i],grid_y[j],grid_z[k]]))
        if abs(cos(dot(g,pom)) - 1.) < 1e-8:    grid_an.append(pom)
        elif abs(cos(dot(g,pom)) + 1.) < 1e-8:  grid_cat.append(pom)
  
  # Setting some king of ionic sizes which will be used in creating the probability distribution
  cat_radius = atomic_scale*min( [ sqrt(dot(grid_cat[0]-grid_cat[j],grid_cat[0]-grid_cat[j])) \
                                   for j in xrange(1,len(grid_cat)) ])
  an_radius  = atomic_scale*min([ sqrt(dot(grid_an[0]-grid_an[j],grid_an[0]-grid_an[j]))\
                                  for j in xrange(1,len(grid_an)) ])

  # Setting the initial prob. distribution
  probability_an  = array([1. for i in range(len(grid_an))])    
  probability_cat = array([1. for i in range(len(grid_cat))])    

  # Which grid points are occupied
  anions_occupy   = []
  cations_occupy  = []

  # distribute anions
  if sum(array(anions))>0:
    for i in range(len(species)):
      if anions[i]==1:
        for j in range(stoichiometry[i]):

          # Placing the first anion
          if len(anions_occupy) == 0:
            # Chose on point on the anion grid
            pom = randint(0,len(grid_an)-1)
            anions_occupy.append(pom)
            
            # Adding the first anion
            result.add_atom = [grid_an[pom],species[i]]
            # Recomputing the new prob distribution for anions
            probability_an  = probability_an - array([prob(x,an_radius,grid_an[pom],result) for x in grid_an]) 
            # Recomputing the new prob distribution for cations due to anions in the cell
            probability_cat = probability_cat - array([prob(x,an_radius,grid_an[pom],result) for x in grid_cat]) 

          # Placing the rest of anions
          else:
            # Setting the probability cutoff
            p_crit = 0.9 
            # Finding points on the anion grid that have the desired probability
            help=[k for k in range(len(grid_an)) if probability_an[k]>p_crit and k not in anions_occupy]
            # In case no points are found reduce the p_crit until you find something
            while len(help)==0:
                p_crit = p_crit-0.05
                help=[k for k in range(len(grid_an)) if probability_an[k]>p_crit and k not in anions_occupy]
            # Repeat the procedure as for the firs anion
            pom = choice(help)
            anions_occupy.append(pom)
            result.add_atom = [grid_an[pom],species[i]]
            probability_an  = probability_an - array([prob(x,an_radius,grid_an[pom],result) for x in grid_an])
            probability_cat = probability_cat - array([prob(x,an_radius,grid_an[pom],result) for x in grid_cat])


  # Distribute cations
  # The same as for anions
  if sum(array(anions))<len(anions):
    for i in range(len(species)):
      if anions[i]==0:
        for j in range(stoichiometry[i]):
          p_crit = 0.9
          help=[k for k in range(len(grid_cat)) if probability_cat[k]>p_crit and k not in cations_occupy]
          while len(help)==0:
            p_crit = p_crit-0.05
            help=[k for k in range(len(grid_cat)) if probability_cat[k]>p_crit and k not in cations_occupy]
          pom = choice(help)
          cations_occupy.append(pom)
          result.add_atom = [grid_cat[pom],species[i]]
          probability_cat = probability_cat - array([prob(x,cat_radius,grid_cat[pom],result) for x in grid_cat])


  # Setting the scale so that the min_distance is fulfilled
  distances=[]
  for atom in result.atoms:
    for n in Neighbors(result,1,atom.pos):
      distances.append(n.distance)

  while result.scale*min(distances) < min_distance:
    result.scale = result.scale + 0.1
    distances=[]
    for atom in result.atoms:
      for n in Neighbors(result,1,atom.pos):
        distances.append(n.distance)

  return result

def fractional_pos(structure, roll=True):
  """ Returns fractional coordinates of atoms in structure. """
  from random import random
  from numpy import dot, array
  from numpy.linalg import inv
  # random translation vector.
  trans = (random(), random(), random()) if roll else (0,0,0)
  result = [dot(inv(structure.cell), atom.pos) + trans for atom in structure]
      
  for pos in result:
    for j in xrange(len(pos)):
      while pos[j] < 0.:  pos[j] += 1. 
      while pos[j] >= 1.: pos[j] -= 1.

  return array(result)

def cut_and_splice(s1, s2, roll=True):
  """ Cut-n-splice GSGO crossover operation

      Mating operation on two parent structures s1 and s2.
      Done on scaled atomic positions (cubic systems) by
      cutting them in half and mixing their upper and 
      lower parts.
  """
  from random import choice, random
  from numpy import dot
  from lada.crystal import Structure

  # swap structures randomly.
  if choice([True, False]): s1, s2 = s2, s1

  # chem. symbols and scaled positions of the two parents
  sc_pos1 = zip([atom.type for atom in s1.atoms],fractional_pos(s1, roll=True))
  sc_pos2 = zip([atom.type for atom in s2.atoms],fractional_pos(s2, roll=True))

  result = Structure()

  result.cell  = s1.cell.copy()
  result.scale = s1.scale

  # choose random positions of split-plane
  xsep = 0.5 - (random() * 0.45 + 0.15)
  # choose direction of split-plane randomly from cell-vectors.
  direction = choice(range(3))

  for type, pos in sc_pos1:
    if pos[direction] >= xsep: result.add_atom = dot(result.cell, pos), type

  for type, pos in sc_pos2:
    if pos[direction] < xsep: result.add_atom = dot(result.cell, pos), type

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
  from ..crystal import Structure

  # swap structures randomly.
  if choice([True, False]): s1, s2 = s2, s1

  # chem. symbols and scaled positions of the two parents
  sc_pos2 = zip([atom.type for atom in s2.atoms],fractional_pos(s2, roll))

  # cell from s1
  result = Structure()
  result.cell  = s1.cell
  result.scale = s1.scale

  # atoms from s2
  for type, pos in sc_pos2:
    result.add_atom = dot(result.cell, pos), type
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
  sc_pos1 = zip([atom.type for atom in s1.atoms], fractional_pos(s1, roll))
  sc_pos2 = zip([atom.type for atom in s2.atoms], fractional_pos(s2, roll))

  result = Structure()
  result.cell  = s1.cell

  for pos, type in chain(sc_pos1, sc_pos2):
    if choice([True, False]):
      result.add_atom = dot(result.cell, type), pos

  result.scale = abs(det(s1.cell)) / abs(det(result.cell)) * s1.scale \
                 * float(len(result.atoms)) / float(len(s1.atoms))

  return result

def jiggle_structure(structure):
  """ Mutates by jiggling structure. """
  from copy import deepcopy
  from numpy import dot, min
  from numpy.linalg import norm
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
  return result


