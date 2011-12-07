""" Tags ground-states. """
__docformat__ = "restructuredtext en"
__all__ = ['convexhull', 'groundstates']

      
def merge_queries(d, filters):
  """ Merge filters into d. """
  if filters is None: return d
  for key, value in filters.iteritems():
    if key in d:
      if isinstance(d[key], dict) and isinstance(value, dict):
        dummy = value.copy(); dummy.update(d[key])
        d[key] = dummy
    else: d[key] = value
  return d

def iter_fere_ternaries(collection="extracted", filters=None, tempname="ladabaseextracteditersystemtempname"):
  """ Loops over FERE ternaries. 
  
      Finds the combination of species which make up a system class in the
      database. To be considered, calculations must be compatible with FERE, eg
      when an enthalpy tag is present in the metadata. Furthermore, multiple
      stoechiometries must be considered.
  """
  from pymongo.code import Code
  from . import Manager
  from .extracted import sort_species
  if isinstance(collection, str):
    ladabase = Manager()
    collection = ladabase.database[collection]

  reduce = Code("function (key, values) { return values[0]; }")
  map = Code( "function() {\n"                                                \
              "  var species = []\n"                                          \
              "  for(specie in this.input.species) species.push(specie)\n"    \
              "  if(species.length == 3) emit(species.toString(), {species: species})\n" \
              "}\n" )
  # map/reduce over objects with FERE enthalpies.
  filters = merge_queries({'metadata.Enthalpy': {'$exists': True}}, filters)
  results = collection.map_reduce(map, reduce, tempname, query=filters)
  # loop over results. This has identified all possible ternary FERE systems.
  for result in results.find(): yield sort_species(result['value']['species'])
  results.drop()


def iter_lowest_energy(species, collection="extracted", filters=None):
  """ Iterates over all lowest-energy structures in a system. 

      Does not check for convex-hull, just for concentration. 
  """ 
  if isinstance(collection, str):
    from . import Manager
    ladabase = Manager()
    collection = ladabase.database[collection]

  # construct query defining the system.
  query = {'metadata.Enthalpy': {'$exists': True}}
  query['metadata.species'] = {'$not': {'$elemMatch': {'$nin': species}}}
  query = merge_queries(query, filters)

  # loops over systems, checking for concentration and enthalpy.
  systems = [ [r['_id'], r['metadata']['Enthalpy'], [r['input']['species'].get(i, 0) for i in species]]\
              for r in collection.find(query) ]
  # construct list with unique concentrations.
  result = []
  while len(systems) > 0:
    result.append(systems.pop())
    Na, Ea, stoecha = sum(result[-1][2]), result[-1][1], result[-1][2]
    for i in xrange(len(systems)-1, -1, -1):
      Nb, Eb, stoechb = sum(systems[i][2]), systems[i][1], systems[i][2]
      if all(a*Nb == b * Na for a, b in zip(stoecha, stoechb)): 
        if Ea > Eb: result[-1] = systems.pop(i)
        else: systems.pop(i)
  
  for r in result: yield collection.find_one({'_id': r[0]})

def half_space_representation(species, collection="extracted", filters=None):
  """ Half-space representation for a given system. 
  
      Returns a 3-tuple consisting of the stoechiometry matrix A,
      the corresponding formation enthalpy vector E, and a identification
      vector listing the corresponding systems in the database.
      The inequalities are set up such that A |Delta||mu| <= E, with
      |Delta||mu| the chemical potential. If the set of inequalities is true,
      the compounds dissociate into their elemental components.


      .. |mu|  unicode:: U+003BC .. GREEK SMALL LETTER MU
      .. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
  """
  from numpy import array
  A, ids = [], []
  for system in iter_lowest_energy(species, collection, filters):
    A.append([system['input']['species'].get(s, 0) for s in species])
    Natoms = 1./float(sum(system['input']['species'].itervalues()))
    A[-1] = [a*Natoms for a  in A[-1]] # goes to concentration.
    A[-1].append(-system['metadata']['Enthalpy'])
    ids.append(system['_id'])
  return array(A), ids

def contour(vertices):
  """ Arranges vertices clockwise around their center. 
  
      It is assumed all vertices are on the same plane.
  """
  from operator import itemgetter
  from numpy import array, sum, arctan2
  if len(vertices) <= 3: return vertices

  center = sum(vertices, axis=0) / float(len(vertices))
  thetas = arctan2( *((vertices-center).T) )
  thetas = sorted(enumerate(thetas), key=itemgetter(1))
  return array([vertices[u[0]] for u in thetas])
 
def generate_fere_summary(filters=3, tempname="ladabaseextracteditersystemtempname"):
  """ Creates FERE ground-states, as well as single-stoechiometry items, in summary database. """
  from datetime import datetime, timedelta
  from numpy import concatenate, identity, zeros, array
  from polyhedron import Hrep
  from . import Manager

  ladabase = Manager()
  extracted = ladabase.database['extracted']
  fere_summary = ladabase.database['fere_summary']

  # Works on last added only
  if isinstance(filters, int): 
    filters = {'metadata.date_added': {'$gte': datetime.now() - timedelta(3)},
               'metadata.Enthalpy': {'$exists': True}}
  elif filters is None: filters = {'metadata.Enthalpy': {'$exists': True}}

  for species in iter_fere_ternaries(extracted, filters, tempname):
  
    # find thermodynamically stable compounds.
    if 'metadata.date_added' in filters: del filters['metadata.date_added']
    A, ids = half_space_representation(species, collection=extracted, filters=filters)
    Nsystems, Nvariables = A.shape[0], A.shape[1]-1
    A = concatenate((A, concatenate((identity(Nvariables), zeros((Nvariables, 1))), axis=1)))
    hrep = Hrep(A[:, :-1], -A[:, -1])

    # find all compounds for this system.
    query = {'metadata.Enthalpy': {'$exists': True},
             'metadata.species': {'$not': {'$elemMatch': {'$nin': species}}}}
    allids = [u['_id'] for u in ladabase.extracted.find(query)]
    document = { '_id': ''.join(species),
                 'species': species,
                 'mus_diagram': { 'formulae': [],
                                  'polyhedra': [],
                                  'extracted_ids': [],
                                  'x': species[0],
                                  'y': species[1]
                            },
                 'extracted_ids': allids,
               }

    # loop over stable compounds only.
    for indices, id in zip(hrep.ininc[:Nsystems], ids):
      if len(indices) == 0: continue
      assert all(array(hrep.is_vertex)[indices] == 1)
      poly = contour(hrep.generators[indices][:, [0, 1]])

      element = extracted.find_one({'_id': id})
      document['mus_diagram']['formulae'].append(element['metadata']['formula'])
      document['mus_diagram']['extracted_ids'].append(element['_id'])
      document['mus_diagram']['polyhedra'].append([u for u in zip(poly[:, 0], poly[:, 1])])

    fere_summary.save(document)
    setme = {'$set': {'tracked.fere_summary': document['mus_diagram']['extracted_ids']}}
    for id in allids: extracted.update({'_id': id}, setme)
