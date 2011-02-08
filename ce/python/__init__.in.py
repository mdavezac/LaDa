""" Cluster Expansion. """
__docformat__ = "restructuredtext en"
from _ce import *
from _fit import *

ce_lattice_type = "@LATTICE_TYPE@"

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
