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
  from numpy import array, concatenate, dot, max
  from polyhedron import Vrep
  individuals = list(individuals)
  # coordinates of individuals
  points = [u.concentrations[:-1].tolist() + [u.energy] for u in individuals]
  # if a convex-hull already exists, adds generators.
  if hasattr(self, 'convexhull'): points.extend(self.convexhull.generators)

  if len(points) == 0: return

  # strip away points at same stiochiometry, keeping only lowest.
  points = strip_points(points)

  # create convex hull from all points.
  vrep = Vrep(points)
  self.convexhull = vrep
  
  # lowers are those half-spaces on the bottom of the energy scale.
  lowers = concatenate((vrep.A, vrep.b[:,None]), axis=1)[vrep.A[:,-1].flat<-0.1]

  # now update individuals
  for indiv in individuals:
    point = array(indiv.concentrations[:-1].tolist() + [indiv.energy])
    indiv.fitness = point[-1] - max(dot(lowers[:,:-2], point[:-1]) - lowers[:,-1])
