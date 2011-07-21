""" Cluster Expansion Module. """
__docformat__ = "restructuredtext en"
from _ce import Cubic, apply_rotation, apply_symmetry, equivalents,\
                ClusterClasses, Clusters, Cluster, create_pairs, \
                find_pis, MLCluster, MLCluster, MLClusterClasses, ce_check
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

def read_jtypes(path, lattice, pairs = 0):
  """ Reads jtype file. """
  from re import compile, M as multline
  from numpy import array
  from ..crystal import lattice_context
  from ..crystal.bravais import bcc, fcc
  from . import MLCluster, MLClusterClasses, _create_clusters

  if lattice == "bcc": lattice = bcc(); lattice.sites[0].type = "A", "B"
  elif lattice == "fcc": lattice = fcc(); lattice.sites[0].type = "A", "B"
  with lattice_context(lattice) as oldlattice:
    if pairs > 0: result = _create_clusters(lattice, 2, pairs) 
    else: result = MLClusterClasses()
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

def cluster_factory(lattice, J0=True, J1=True, **mb):
  """ Returns class of cluster classes. 

      :Parameters:
        lattice : `lada.crystal.Lattice`
          The backbone lattice on which the Ising model is created.
        J0 : boolean
          If True, adds zero order term.
        J1 : boolean
          If True, adds on-site terms.
        **mb: 
          All other many body terms, where the keys should be "B2", "B3", "Bn",
          where n is the order of the many-body interaction, and the values is
          the number of shells to look for.
  """
  from re import compile, match
  from ._ce import _create_clusters
  from ..crystal import which_site

  # checks that keywords are well formed.
  key_regex = compile("B(\d+)")
  for key in mb.keys(): 
    a_re = match(key_regex, key)
    assert a_re != None, "Keyword %s is not of the form B(\d+)" % (key)
    assert a_re.end() == len(key), "Keyword %s is not of the form B(\d+)" % (key)
    assert int(a_re.group(1)) > 1, "Cannot understand keyword %s" % (key)
    assert mb[key] > 0, "Cannot understand input %s=%i" % (key, mb[key])

  # sanity check.
  assert len(mb) > 0 or J0 or J1, ValueError("No clusters to create.")

  # computes equivalent sites.
  equiv_sites = set( i for i in range(len(lattice.sites)) ) 
  for i, site in enumerate(lattice.sites):
    if i not in equiv_sites: continue
    for op in lattice.space_group:
      j = which_site( op(site.pos), lattice )
      if j != i and j in requiv_sites: equiv_sites.remove(j)

  result = None
  # now creates multi-lattice clusters, along with index bookkeeping.
  if J0: result = _create_clusters(lattice, nth_shell=0, order=0, site=0) # J0
  # creates J1 for these sites.
  if J1: 
    for site in equiv_sites:
      dummy = _create_clusters(lattice, nth_shell=0, order=1, site=site)
      if result == None: result = dummy
      else: result.extend(dummy)
  # creates many bodies.
  for site in equiv_sites:
    for key in sorted(mb.keys(), cmp):
      regex = match(key_regex, key)
      order = int(regex.group(1))
      shell = mb[key]
      dummy = _create_clusters(lattice, nth_shell=shell, order=order, site=site)
      if result == None: result = dummy
      else: result.extend(dummy)
  return result
