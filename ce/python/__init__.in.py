""" Cluster Expansion Module. """
__docformat__ = "restructuredtext en"
from _ce import Cubic, apply_rotation, apply_symmetry, equivalents,\
                ClusterClasses, Clusters, Cluster, create_pairs, create_clusters, \
                find_pis, MLCluster, \
                MLCluster as _MLClusters, MLClusterClasses as _MLClusterClasses
from _fit import *

def _lattdeco(method):
  """ Switches to self-owned lattice for duration of call. """
  def wrapped(self, *args, **kwargs):
    from ..crystal import lattice_context
    with lattice_context(self.lattice) as oldlattice:
      return method(self, *args, **kwargs)
  wrapped.__name__ = method.__name__
  wrapped.__doc__ = method.__doc__
  wrapped.__module__ = method.__module__
  if hasattr(method, '__dict__'): wrapped.__dict__.update(method.__dict__)
  return wrapped

ce_lattice_type = "@LATTICE_TYPE@"

class MLClusters(_MLClusters):
  """ Array of equivalent clusters.
  
      :note: This implementation ties the cluster to a back-bone lattice. The
      global lattice is switched to that lattice for the duration of calls.
  """
  def __init__(self, lattice, cluster):
    """ Initializes a class of equivalent clusters.

        :Parameters:
          lattice : `Lattice`
            Lattice for which to create cluster
          cluster : `MLCluster`
            Cluster from which to create a class of equivalent clusters.
    """
    from ..crystal import lattice_context
    self._lattice = lattice
    """ Backbone lattice on which clusters reside. """
    with lattice_context(self.lattice) as dummy: 
      super(MLClusters, self).__init__(cluster)
  @property
  def lattice(self):
    """ Backbone lattice on which cluster resides. """
    return self._lattice

  extend   = _lattdeco(_MLClusterClasses.extend)

   
class MLClusterClasses(_MLClusterClasses):
  """ Array of `MLClusters`. 
 
      :note: This implementation ties the cluster to a back-bone lattice. The
      global lattice is switched to that lattice for the duration of calls.
  """

  def __init__(self, clusters = None, lattice = None):
    from . import MLCluster
    from ..crystal import lattice_context
    assert lattice != None or input != None, \
           ValueError('One of clusters or lattice must given on input.')
    if lattice != None: self._lattice = lattice
    if clusters == None: return
    if lattice == None:
      clusters = [u for u in clusters]
      # finds lattice
      for u in clusters: 
        if hasattr(u, 'lattice'):
          self._lattice = u.lattice 
          break
    assert self.__dict__.get('_lattice', None) != None,\
           ValueError("Could not determine lattice from input.")
    with lattice_context(self.lattice) as previous:
      super(MLClusterClasses, self).__init__([MLCluster(u) for u in clusters])

  @property
  def lattice(self):
    """ Backbone lattice on which classes of clusters reside. """
    return self._lattice

  append   = _lattdeco(_MLClusterClasses.append)
  extend   = _lattdeco(_MLClusterClasses.extend)
  pis      = _lattdeco(_MLClusterClasses.pis)
  __call__ = _lattdeco(_MLClusterClasses.__call__)




def read_jtypes(path, lattice, pairs = 0):
  """ Reads jtype file. """
  from re import compile, M as multline
  from numpy import array
  from ..crystal import lattice_context
  from ..crystal.bravais import bcc, fcc
  from . import MLCluster, MLClusters, MLClusterClasses, create_clusters

  if lattice == "bcc": lattice = bcc(); lattice.sites[0].type = "A", "B"
  elif lattice == "fcc": lattice = fcc(); lattice.sites[0].type = "A", "B"
  with lattice_context(lattice) as oldlattice:
    if pairs > 0: result = create_clusters(lattice, 2, pairs) 
    else: result = MLClusterClasses(lattice)
    J0 = compile(r"^\s*J0(?:\s|#)*$\s*0\s+0(?:\s|#)*$", multline)
    J1 = compile(r"^\s*J1(?:\s|#)*$\s*1\s+1(?:\s|#)*$\s*0\s+0\s+0(?:\s|#)*$", multline)
    mb = compile(r"^\s*(\dB\d-\d+)(?:\s|#)*$"\
                 r"\s*\d+\s+\d+(?:\s|#)*$"\
                 r"(?:\s*\d+\s+\d+\s+\d+(?:\s|#)*)*$", multline)
  
    with open(path, "r") as file: text = file.read()
    if J0.search(text) != None: result.append()
    if J1.search(text) != None: result.append(MLCluster())
    for manybody in mb.finditer(text):
      cluster = MLCluster()
      cluster.origin.site = 0
      cluster.origin.pos = array([0,0,0], dtype="float64")
      for vector in manybody.group(0).split('\n')[3:]:
        cluster.append([float(u) * 0.5 for u in vector.split()], 0)
      result.append(cluster)

    return result


def read_mbce_structure(path):
  """ Read MBCE structure. """
  from numpy import array, zeros
  from ..crystal import Structure
  structure = Structure()
  with open(path, "r") as file: lines = file.readlines()
  structure.name = lines[0].rstrip().lstrip()
  cell = zeros((3,3), dtype="float64")
  for i in range(3): cell[:,i] = array(lines[i+2].split(), dtype="float64")
  structure.cell = cell
  for line in lines[5:]:
    if len(line.split()) < 4: break
    if line.split()[0] not in ['1', '2']: continue 
    structure.add_atom = array(line.split()[1:4], dtype="float64"),\
                         "A" if line.split()[0] == "1" else "B"
  structure.scale  = 1e0
  structure.energy = 0e0
  return structure

def read_mbce_structures(directory):
  """ Reads a list of MBCE structures. """
  from os.path import exists, join
  assert exists(join(directory, "LDAs.dat")), RuntimeError("LDAs.dat file could not be found.")
  result = []
  with open(join(directory, "LDAs.dat"), 'r') as file:
    for line in file:
      if line.find('#') != -1: line = line[:line.find('#')]
      if len(line.split()) < 2: continue
      result.append( read_mbce_structure(join(directory, line.split()[0])) )
      result[-1].energy = float(line.split()[1])
  return result
