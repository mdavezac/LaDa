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

import numpy
from lada.crystal import A2BX4
from lada.opt import read_input
from lada.ce import create_clusters

# input
shell = 5

# create spinel lattice.
lattice = A2BX4.b5()
# occupations.
for site in lattice.sites:
  if "X" in site.type: continue
  site.type = ["A", "B"]
# recomputes space group for safety.
lattice.find_space_group()
lattice.set_as_crystal_lattice()

# list of inequivalent sites, keeping only those without "X"
ineqs = []
for i in  inequivalent_sites(lattice):
  if len(lattice.sites[i].type) > 1: ineqs.append(i)

# now creates multi-lattice clusters, along with index bookkeeping.
mlclasses = create_clusters(lattice, nth_shell=0, order=0, site=0) # J0
print len(mlclasses)
# creates J1 for these sites.
for site in ineqs:
  mlclasses.extend(create_clusters(lattice, nth_shell=0, order=1, site=site))
# create pairs clusters.
for site in ineqs:
  mlclasses = create_clusters(lattice, nth_shell=shell, order=2, site=site)
print len(mlclasses)

