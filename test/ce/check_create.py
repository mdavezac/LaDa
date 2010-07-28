def inequivalent_sites(lattice):
  """ Returns a list containing only one site index per inequivalent sub-lattice. """
  from lada.crystal import which_site

  result = set( i for i in range(len(lattice.sites)) ) 
  for i, site in enumerate(lattice.sites):
    if i not in result: continue
    for op in lattice.space_group:
      j = which_site( op(site.pos), lattice )
      if j != i and j in result: result.remove(j)
  return result

from numpy import array, zeros
from numpy.linalg import norm
import pyublas 
from random import choice
from lada.crystal import A2BX4, fill_structure
from lada.opt import read_input
from lada.ce import create_clusters
from math import exp

# input
shell = 35

# create spinel lattice.
lattice = A2BX4.b5()
# occupations.
for site in lattice.sites:
  if "X" in site.type: continue
  site.type = ["A", "B"]
# recomputes space group for safety.
lattice.find_space_group()
lattice.set_as_crystal_lattice()

# creates random structure.
structure = fill_structure(A2BX4.b5I().cell, A2BX4.b5I())
normal = fill_structure(A2BX4.b5().cell, A2BX4.b5())
# for atom in structure.atoms:
#   atom.type = choice(lattice.sites[atom.site].type)


# list of inequivalent sites, keeping only those without "X"
ineqs = []
for i in  inequivalent_sites(lattice):
  if len(lattice.sites[i].type) > 1: ineqs.append(i)


# now creates multi-lattice clusters, along with index bookkeeping.
# first J0
clusters = create_clusters(lattice, nth_shell=0, order=0, site=0)
# then J1 for these sites.
for i in ineqs:
  clusters.extend(create_clusters(lattice, nth_shell=0, order=1, site=i))
# then pairs figures.
clusters = None
for i in ineqs:
  if clusters == None:
    clusters = create_clusters(lattice, nth_shell=shell, order=2, site=i)
  else: 
    clusters.extend(create_clusters(lattice, nth_shell=shell, order=2, site=i))
# sets interaction energies (eci)
# J0
clusters[0].eci = 0e0
# J1
for i in range(1, len(ineqs)+1): clusters[i].eci = 0e0
# J2
# clusters -> list of pair figures.
# clusters[i] -> a list of symmetrically *equivalent* figures.
# clusters[i][0] -> one single figure.
# clusters[i][0].origin -> origin of the figure: eg first spin.
# clusters[i][0][j] -> vector from origin to other spins:
#    eg clusters[i][0].origin.pos + clusters[i][0][0].pos = position of second spin
for i in range(len(ineqs)+1, len(clusters)): 
# print norm(clusters[i][0][0].pos)
  clusters[i].eci = exp(-5e0*norm(clusters[i][0][0].pos))

# now computes pis.
pis = clusters.pis(structure) - clusters.pis(normal)
# print pis
# now compute energies.
print shell, clusters(structure) - clusters(normal)
