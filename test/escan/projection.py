import pickle
import matplotlib.pyplot as plt 
from lada.escan import ExtractBS

def compute_projection(extract, alpha = 1e0, **kwargs):
  """ Computes projected densities around each atom for a calculation. """
  from os.path import join
  from pickle import dump
  from numpy import zeros, dot, array
  from numpy.linalg import norm
  from lada import periodic_table as table
  from lada.crystal import to_voronoi
  species = set([u.type for u in extract.vff.structure.atoms])
  species = sorted(list(species))
  for key in species:
    if (key not in kwargs) and (key not in table): continue
    radius = kwargs.get(key, getattr(table, key))
    radius *= radius
    sigma = -alpha / radius

    # create projection operator.
    proj = zeros(extract.rvectors.shape)
    cell = extract.vff.structure.cell
    for atom in extract.vff.structure.atoms:
      if atom.type != key: continue
      centered = to_voronoi(extract.rvectors, atom.pos, cell)
      proj += exp( dot(centered, centered.T) * sigma )

    result[key] = [(w.eigenvalue, w.expectation_value(proj)) for w in extractor.rwfns]
  
  n = norm(array(result.values()))
  for key in result.keys(): result[key] /= n
  with open(join(extract.directory, "PROJECT_BS")) as file: dump(result, file)

def read_projections(extract):
  from os.path import join
  from pickle import load
  with open(join(extract.directory, "PROJECT_BS")) as file: return load(file)

import pickle
from sys import exit
from os.path import join
from numpy import matrix, array
from numpy.linalg import norm
from boost.mpi import world
from lada.opt import read_input
from lada.escan import Escan, soH, band_structure
from lada.vff import Vff
from lada.crystal import fill_structure, sort_layers, FreezeCell, nb_valence_states  

# reads input file.
global_dict={"Vff": Vff, "Escan": Escan, "nb_valence_states": nb_valence_states, "soH": soH}
input = read_input("input.py", global_dict=global_dict)

# creating unrelaxed structure.
structure = input.vff.lattice.to_structure()
structure.atoms[0].type = "Si"
structure.atoms[1].type = "Ge"
structure.scale = 5.65

# some kpoints + associated name
X = array( [1,0,0], dtype="float64" )
G = array( [0,0,0], dtype="float64" )
L = array( [0.5,0.5,0.5], dtype="float64" )
W = array( [0, 0.5,1], dtype="float64" )

# Each job is performed for a given kpoint (first argument), at a given
# reference energy (third argument). Results are stored in a specific directory
# (second arguement). The expected eigenvalues are given in the fourth argument.
kpoints = [ (X, G), (G, L) ]
density = 20 / min( norm(X), norm(L), norm(W) )

result = band_structure( input.escan, structure, kpoints, density, 
                         outdir = "results",
                         eref   = None, 
                         nbstates = nb_valence_states(structure) + 4,
                         pools = 4)
  
if world.rank == 0:
  with open(join("results", "pickle"), "w") as file:
    pickle.dump(result, file) 

