""" Tags ground-states. """
__docformat__ = "restructuredtext en"
__all__ = ['convexhull', 'groundstates']

def add_elements(collection='extracted'):
  """ Adds elemental energies. """
  from itertools import chain
  from quantities import eV
  from . import Manager
  mus_final = { 'Ag': -0.82700958541595615*eV,
                'Au': -2.2303410086960551*eV,
                'Ba': -1.3944992462870172*eV,
                'Be': -3.3972092621754264*eV,
                'Bi': -4.3853003286558812*eV,
                'Cl': -1.6262437135301639*eV,
                'Co': -4.7543486260270402*eV,
                'Cr': -7.2224146752384204*eV,
                'Cu': -1.9725806522979044*eV,
                'F' : -1.7037867766570287*eV,
                'Fe': -6.1521343161090325*eV,
                'Ge': -4.137439286830797*eV,
                'Hf': -7.397695761161847*eV,
                'Hg': -0.12361566177444684*eV,
                'Ir': -5.964577394407752*eV,
                'K' : -0.80499202755075006*eV,
                'La': -3.6642174822805287*eV,
                'Li': -1.6529591953887741*eV,
                'Mn': -6.9965778258511993*eV,
                'Na': -1.0640326227725869*eV,
                'Nb': -6.6867516375690608*eV,
                'Ni': -3.5687859474688026*eV,
                'Pd': -3.1174044624888873*eV,
                'Pt': -3.9527597082085424*eV,
                'Rb': -0.6750560483522855*eV,
                'Rh': -4.7622899695820369*eV,
                'Sb': -4.2862260747305099*eV,
                'Sc': -4.6302422200922519*eV,
                'Si': -4.9927748122726356*eV,
                'Sn': -3.7894939351245469*eV,
                'Sr': -1.1674559193419329*eV,
                'Ta': -8.8184831379805324*eV,
                'Te': -3.2503408197224912*eV,
                'Ti': -5.5167842601434147*eV,
                'V' : -6.4219725884764864*eV,
                'Y' : -4.812621315561298*eV,
                'Zr': -5.8747056261113126*eV}

  mus_fix = {'O' :-4.76*eV,
             'S' :-4.00*eV,
             'Se':-3.55*eV,
             'N' :-8.51*eV,
             'P' :-5.64*eV,
             'As':-5.06*eV,
             'Mg':-0.99*eV,
             'Ca':-1.64*eV,
             'Zn':-0.84*eV,
             'Cd':-0.56*eV,
             'Ga':-2.37*eV,
             'Al':-3.02*eV,
             'In':-2.31*eV}

  ladabase = Manager()
  collection = ladabase.database[collection]
  for key, value in chain(mus_final.iteritems(), mus_fix.iteritems()):
    indatabase = [k for k in collection.find({'elemental': key})]
    r = {'total_energy': value.magnitude.tolist(), 'elemental': key}
    if len(indatabase) == 1: r['_id'] = indatabase[0]['_id']
    elif len(indatabase) > 1: 
      raise RuntimeError('found more than one items corresponding to the element {0}.'.format(key))
    collection.save(r)
      
def merge_queries(d, filters):
  """ Merge filters into d. """
  if filters == None: return d
  for key, value in filters.iteritems():
    if key in d:
      if isinstance(d[key], dict) and isinstance(value, dict):
        dummy = value.copy(); dummy.update(d[key])
        d[key] = dummy
    else: d[key] = value
  return d

def convexhull(species, collection="extracted", filters=None, withargs=False, tolerance=0.003, verbose=False):
  """ Determines convex-hull of given specie set. """
  from numpy import multiply, sum, array
  from ..hull import Hull3d
  from . import Manager
  assert len(species) == 3
  
  if isinstance(collection, str):
    ladabase = Manager()
    collection = ladabase.database[collection]

  # list of all systems with elements of interest.
  query = merge_queries({'species': {'$all': species}}, filters)
  fields = ['total_energy', 'species', 'ions_per_specie', 'system']
  systems = [k for k in collection.find(query, fields=fields)]

  # correctly ordered list of elemental energies.
  elementals = [collection.find_one({'elemental': k})['total_energy'] for k in species]
  # loop over systems an add a concentration.
  def point3d(group, stoechiometry, energy):
    """ Computes concentration from species and stoeckiometry. """
    stoechiometry = array([(0 if u not in group else float(stoechiometry[group.index(u)])) for u in species])
    energy -= sum(multiply(elementals, stoechiometry))
    N = float(sum(stoechiometry))
    return [stoechiometry[0] / N, stoechiometry[1] / N, energy / N]

  points = [ point3d(system['species'], system['ions_per_specie'], system['total_energy']) \
             for system in systems ]
  if verbose: print "  Found {0} points, excluding elementals.".format(len(points))
  points.append([1, 0, 0])
  points.append([0, 1, 0])
  points.append([0, 0, 0])

  # creates and returns convex hull.
  hull = Hull3d(points, tolerance=tolerance**2)
  if not withargs: return hull
  for system, point in zip(systems, points):
    system['distance'] = hull.distance(point)
    system['deltaH'] = point[2]
  return hull, systems


def groundstates(collection="extracted", field="groundstate", tolerance=0.003, filters=None):
  """ Trolls database to find all groundstates. """
  from time import time
  from numpy import abs
  from . import Manager

  ladabase = Manager()
  collection = ladabase.database[collection]

  # set of elemental energies.
  elementals = set([k['elemental'] for k in collection.find({'elemental' : {'$exists': 1}})])
  
  # creates list of unique groups with three, if the elemental values are known.
  query = merge_queries({'species': {'$size': 3}, 'total_energy': {'$exists': 1}}, filters)
  threeiter = collection.find(query, fields=['species'])
  threes = []
  for group in threeiter:
    group = set(group['species'])
    # check it is not already counted.
    if group in threes: continue
    # check elemental values are known.
    if group.issubset(elementals): threes.append(group)
    else: print group

  # now loop over 3d convex-hulls.
  for group in threes:
    group = list(group)
    print "Working on {0}.".format(group)
    timing = time()
    hull, systems = convexhull( group, collection=collection, filters=filters,
                                withargs=True, tolerance=tolerance, verbose=True )
    timing = time() - timing
    hour = int(float(timing/3600e0))
    minute = int(float((timing - hour*3600)/60e0))
    second = (timing - hour*3600-minute*60)
    print "  created-convex hull in {0:2>}:{0:2>}:{0:12.8e}.".format(hour, minute, second)
    for system in systems:
      if abs(system['distance']) <= tolerance: 
        collection.update( {'_id': system['_id']},
                           {'$set': {field: 1, 'distance': system['distance'], 'deltaH': system['deltaH']}})
      else:
        collection.update( {'_id': system['_id']},
                           { '$unset': {field: 1},
                             '$set': {'distance': system['distance'], 'deltaH': system['deltaH']} } )
    
  # no twos yet.
  # twoiter = collection.find( { 'species': {'$size': 2},
  #                              'total_energy': {'$exists': 1} },
  #                            fields=['species'] )
  # twos = []
  # for group in twoiter:
  #   group = set(group['species'])
  #   # check it is not already counted.
  #   if group in twos: continue
  #   # check elemental values are known.
  #   if not group.issubset(elementals): continue
  #   # check that is not already included in threes.
  #   ok = True
  #   for three in threes: 
  #     if group.issubset(three): ok = False; break
  #   if ok: twos.append(group)

    
