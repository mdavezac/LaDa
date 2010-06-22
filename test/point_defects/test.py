# will be params:
from sys import exit
from os.path import join
from lada.opt import read_input
from lada.vasp import Vasp, Specie, specie
from lada.vasp.incar import Standard, NElect
from lada.crystal import fill_structure, Neighbors, point_defects as ptd
from lada.crystal.
import jobs

# names we need to create input.
input_dict = {
               "Specie": Specie, "Vasp": Vasp, "Standard":Standard,
               "U": specie.U, "nlep": specie.nlep
             }
# reads input.
input = read_input("input.py", input_dict)
# sets as default lattice.
input.lattice.set_as_crystal_lattice()
# changes species in lattice.
for site in input.lattice.sites:
  if   "A" in site.type: site.type[0] = input.species[0].symbol
  elif "B" in site.type: site.type[0] = input.species[1].symbol
  elif "X" in site.type: site.type[0] = input.species[2].symbol

# creates super-structure
supercell = fill_structure(input.supercell)

# vacancy
for symbol, specie in input.species.items(): 
  for structure, vacancy in ptd.vacancy(supercell, input.lattice, symbol):
    vacdir = join(outdir, "vacancy_" + symbol)

    # max/min oxidation state
    max_oxidation = specie.oxydation if hasattr(specie, "oxydation") else 0
    # sign of oxidation state
    sign_oxidation = 1 if max_oxidation > 0 else -1
    # loop over oxidation states.
    for oxidation in range(0, max_oxidation + sign_oxidation, sign_oxidation): 
      oxdir = join(vacdir, str(i) + "-" if oxidation > 0 else "+") if max_oxidation != 0 else vacdir

      # now finds first neighbors. 12 is the highest coordination number, so
      # this should include the first shell.
      neighbors = [n for n in Neighbors(structure, 12, vacancy.pos)]
      # only take the first shell.
      neighbors = [n for n in neighbors if n.distance < neighbors[0].distance + 1e-12]
      # reduce to those which are magnetic.
      neighbors = [n for n in neighbors if input.vasp.species[ structure.atoms[n.index] ].magnetic]
      # reduces to sorted list of indices in structure.
      indices = sorted([n.index for n in neighbors])
      # creates a list of species.
      species = [input.vasp.species[ structure.atoms[n] ] for n in neighbors]

      if len(neighbors) == 0: # no magnetic neighbors.
        jobs.current[oxdir].job["structure"] = structure
        jobs.current[oxdir].job["nelect"] = NElect(-oxidation)
        jobs.current[oxdir].job["vasp"] = input.vasp
      else: # Creates low spin and high-spin *ferro* configurations.

        # First, find low-spin value of each magnetic ion. Assumes that
        # magnetic species are s^2 d^n p^0!
        lowspin = [s.valence % 2  + float(oxidation) / float(len(species)) for s in species]
        magmom, last_index = "", 0
        for index, valence in zip(indices, lowspin):
          if index - last_index == 1:   magmom += "%f " % (valence)
          elif index - last_index == 2: magmom += "0 %f " % (valence)
          else:                         magmom += "%i*0 %f" % (index-last_index-1, valence)
          last_index = index
        
        lowdir = join(oxdir, "lowspin")
        jobs.current[lowdir].job["structure"] = structure
        jobs.current[lowdir].job["nelect"] = NElect(-oxidation)
        jobs.current[lowdir].job["magmom"] = magmom
        jobs.current[lowdir].job["vasp"] = input.vasp
        

        # Then, find low-spin value of each magnetic ion. Assumes that
        # magnetic species are s^2 d^n p^0! Moment may be too high here.
        mag = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0]
        highspin = [mag[s.valence-2] for s in species]
        magmom, last_index = "", 0
        for index, valence in zip(indices, highspin):
          if index - last_index == 1:   magmom += "%f " % (valence)
          elif index - last_index == 2: magmom += "0 %f " % (valence)
          else:                         magmom += "%i*0 %f" % (index-last_index-1, valence)
          last_index = index
        
        highdir = join(oxdir, "highpin")
        jobs.current[highdir].job["structure"] = structure
        jobs.current[highdir].job["nelect"] = NElect(-oxidation)
        jobs.current[highdir].job["magmom"] = magmom
        jobs.current[highdir].job["vasp"] = input.vasp
        



      



    # case with no oxidation states.
    if not hasattr(vac, "oxidation"): 
      jobs.current[vacdir].job["structure"] = structure
      jobs.current[vacdir].job["vasp"] = input.vasp
    # case with oxidation states.
    else:
      for i in range(0, vac.oxidation, 1 if vac.oxidation > 0 else -1): 
        directory = join(vacdir, str(i) + "-" if vac.oxidation > 0 else "+")
        jobs.current[directory].job["structure"] = structure
        jobs.current[directory].job["nelect"] = NElect(-i)
        jobs.current[directory].job["vasp"] = input.vasp

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


