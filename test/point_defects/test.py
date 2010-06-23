# will be params:
from sys import exit
from os.path import join
from lada.opt import read_input
from lada.vasp import Vasp, Specie, specie, files
from lada.crystal import fill_structure, Neighbors, point_defects as ptd
from lada import jobs

# names we need to create input.
input_dict = { "Vasp": Vasp, "U": specie.U, "nlep": specie.nlep }
# reads input.
input = read_input("input.py", input_dict)
# sets as default lattice.
input.lattice.set_as_crystal_lattice()

# creates super-structure
supercell = fill_structure(input.supercell)

# loop over specie vacancy
for symbol, specie in input.vasp.species.items(): 
  # loop over types of vacancies for given specie (eg symmetrically inequivalent sites).
  for structure, vacancy in ptd.vacancy(supercell, input.lattice, symbol):
    vacdir = "vacancy_" + symbol
    # max/min oxidation state
    max_oxidation = specie.oxidation if hasattr(specie, "oxidation") else 0
    # sign of oxidation state
    sign_oxidation = 1 if max_oxidation > 0 else -1
    # loop over oxidation states.
    for oxidation in range(0, max_oxidation + sign_oxidation, sign_oxidation): 
      # directory
      if max_oxidation == 0: oxdir = vacdir
      elif oxidation == 0:   oxdir = join(vacdir, "neutral")
      elif oxidation > 0:    oxdir = join(vacdir, "+" + str(oxidation) )
      elif oxidation < 0:    oxdir = join(vacdir, str(oxidation) )

      # now finds first neighbors. 12 is the highest coordination number, so
      # this should include the first shell.
      neighbors = [n for n in Neighbors(structure, 12, vacancy.pos)]
      # only take the first shell and keep indices (to atom in structure) only.
      neighbors = [n.index for n in neighbors if n.distance < neighbors[0].distance + 1e-1]
      # reduce to those which are magnetic.
      neighbors = [n for n in neighbors if input.vasp.species[ structure.atoms[n].type ].magnetic]

      if len(neighbors) == 0: # no magnetic neighbors.
        jobs.current[oxdir].job["structure"] = structure
        jobs.current[oxdir].job["nelect"] = -oxidation
        jobs.current[oxdir].job["vasp"] = input.vasp
      else: # has magnetic neighbors
        # Creates low spin and high-spin *ferro* configurations.
        # The current implementation assumes that magnetic species are s^2 d^n p^0!
        # loops over high and low spin configurations.
        for spin in ["low", "high"]:
          spindir = join(oxdir, spin + "-spin")
          jobs.current[spindir].job["structure"] = structure
          jobs.current[spindir].job["nelect"] = -oxidation
          jobs.current[spindir].job["magmom"] = spin, neighbors
          jobs.current[spindir].job["vasp"] = input.vasp

# for job, name in jobs.walk_through(outdir="results"):
#   print name
#   job.compute(norun=True, repat=files.input, outdir=name)
# jobs.current.update(jobs.load())
# for job in jobs.walk_through():
#   print job.job["outdir"]
#   if "nelect" in job.job: print job.job["nelect"]

# # substitutions.
# subs = [(input.species[0].symbol, input.species[1].symbol),\
#         (input.species[1].symbol, input.species[0].symbol)]
# for A, B in subs:
#  for structure in substitution(supercell, input.lattice, A, B):
#    directory = join(outdir, "%s_on_%s" % (A, B))
#    print directory, len([0 for atom in structure.atoms if atom.type == B])


