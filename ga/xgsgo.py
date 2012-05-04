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

def rand_struct_gen( species, stoichiometry, anions=None, angle_range=None,
                     min_distance=1.7, atomic_scale=2., cubic=False):
  """ Generates random crystal structure.
      
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
      :param angle_range: 
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
  from random import randint, shuffle, choice
  from numpy import pi, sqrt, cos, dot, array, arange, ones,      \
                    argwhere, any, exp, sum, logical_not, argmin, \
                    identity, outer, argmax
  from numpy.linalg import inv, norm
  from itertools import product, chain
  from lada.crystal import Structure, neighbors, coordination_shells

  if anions == None: anions = [True]*len(species)
  anions = array(anions, dtype="bool")
  species = array(species)
  stoichiometry = array(stoichiometry, dtype="int")

  minvol = 4.0/3.0*atomic_scale**3 *sum(stoichiometry)
  result = Structure(random_cell(angle_range, cubic))
  result.scale = (2*minvol/result.volume)**(1.0/3.0)

  print result.volume, sum(stoichiometry)
  # Setting two grids, one for the anions and one for the cations
  # This is done creating two grids offset from one another.
  # choice between 4 and 7 is to have large enough number of points, but not
  # too large grids
  N = max(7, (3.0*sum(stoichiometry))**(1./3.))
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
  # now adds each atom in turn. Species do not matter yet.
  for i in xrange(sum(stoichiometry)):
    # create criteria for placing atoms.
    maxi, mini = max(probability), min(probability)
    criteria = 0.90 * (maxi - mini) + mini
    # Find points on the anion grid that have the desired probability
    empty_space=argwhere(probability >= criteria) 
    print len(empty_space)
    # choose position randomly
    assert len(empty_space)
    pom = choice(empty_space.flatten())
    # add atom to structure
    result.add_atom(*grid[pom].flatten(), type='anion')
    # make neighbors less probable
    distances = (grid- grid[pom])[:,None, :] - periodic_images[None, :, :] 
    probability -= sum(exp(-sum(distances**2, axis=2)/atomic_scale**2), axis=1)

  # Set the scale so that the min_distance is fulfilled
  mindist = [ (vector, distance) for atom in result \
              for neighbor, vector, distance in neighbors(result, 1, atom) ]
  mindist = mindist[argmin([d[1] for d in mindist])]
  if mindist[1] * result.scale < min_distance:
    result.scale *= min_distance / result.scale / mindist[1] 

  # now plop down cations 
  cations = logical_not(anions)
  if all(cations):  # no thinking involved in this case.
    for a in result: a.type = 'cation'
  # try and make sure there are as few possible cation-cation and anion-anion
  # neighbors.
  elif any(cations):
    # shell of each atom
    shells = [ coordination_shells(result, 1, atom, 0.4/result.scale) \
               for atom in result ]
    # loop over sum of cations
    for i in xrange(sum(stoichiometry[cations])):
      environment = []
      for shell, atom in zip(shells, result):
        x = [(1.0 if atom.type == n.type else 0.) for n, v, d in shell]
        environment.append(sum(x))
      environment = array(enviromnent)
      result[argmax(environment)].type = 'cation'

  # now distribute cations.
  allcations = ( [s]*n for s, n\
                 in zip(species[cations], stoichiometry[cations]))
  allcations = list(chain(*allcations))
  shuffle(allcations)
  for atom in result:
    if result.type == 'cation': result.type = allcations.pop()
  # now distribute anions.
  allanions = ( [s]*n for s, n\
                in zip(species[anions], stoichiometry[anions]))
  allanions = list(chain(*allanions))
  shuffle(allanions)
  for atom in result:
    if result.type == 'anion': result.type = allanions.pop()

  return result

