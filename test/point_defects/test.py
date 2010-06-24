# will be params:
import cPickle
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
  vacjobs = jobs.current / ("vacancy_" + symbol)
  # loop over types of vacancies for given specie (eg symmetrically inequivalent sites).
  for i, (structure, vacancy) in enumerate(ptd.vacancy(supercell, input.lattice, symbol)):
    vacjob = vacjobs / ("site_" + str(i))

    # loop over oxidation states.
    for oxidation, name in ptd.oxidation(specie):
      oxjob = vacjob / name

      # now finds first neighbors. 12 is the highest coordination number, so
      # this should include the first shell.
      neighbors = [n for n in Neighbors(structure, 12, vacancy.pos)]
      # only take the first shell and keep indices (to atom in structure) only.
      neighbors = [n.index for n in neighbors if n.distance < neighbors[0].distance + 1e-1]
      # reduce to those which are magnetic.
      neighbors = [n for n in neighbors if input.vasp.species[ structure.atoms[n].type ].magnetic]

      if len(neighbors) == 0: # no magnetic neighbors.
        oxjob.add_param = "structure", structure
        oxjob.add_param = "nelect",    -oxidation
        oxjob.vasp = input.vasp
        # Note: structure and nelect can now be accessed as
        # oxjob.structure and oxjob.nelect
      else: # has magnetic neighbors
        # Creates low spin and high-spin *ferro* configurations.
        # The current implementation assumes that magnetic species are s^2 d^n p^0!
        # loops over high and low spin configurations.
        for spin in ["low", "high"]:
          spinjob = oxjob / (spin + "-spin")
          spinjob.add_param = "structure",   structure
          spinjob.add_param = "nelect",      -oxidation
          spinjob.add_param = "magmom",      (spin, neighbors)
          spinjob.vasp = input.vasp
          # Note: structure, magmom, and nelect can now be accessed as
          # spinjob.structure, spinjob.magmom, and spinjob.nelect

# loop substitutions.
for A, B in [("Rh", "Zn"), ("Zn", "Rh") ]:
  # loop over inequivalent substitution sites.
  subjobs = jobs.current / (A + "_on_" + B)
  for i, (structure, substitution) in enumerate(ptd.substitution(supercell, input.lattice, A, B)):
    subjob = subjobs / ("site_" + str(i))

    # loop over oxidations.
    for oxidation, name in ptd.oxidation(input.vasp.species[A], input.vasp.species[B]):
      oxjob = subjob / name

      # now finds first neighbors. 12 is the highest coordination number, so
      # this should include the first shell.
      neighbors = [n for n in Neighbors(structure, 12, vacancy.pos)]
      # only take the first shell and keep indices (to atom in structure) only.
      neighbors = [n.index for n in neighbors if n.distance < neighbors[0].distance + 1e-1]
      # adds origin.
      neighbors.insert(substitution.index, 0)
      # reduce to those which are magnetic.
      neighbors = [n for n in neighbors if input.vasp.species[ structure.atoms[n].type ].magnetic]

      if len(neighbors) == 0: # no magnetic neighbors.
        oxjob.jobparams["structure"] = structure
        oxjob.jobparams["nelect"]    = -oxidation
        oxjob.vasp = input.vasp
      else: # has magnetic neighbors
        # Creates low spin and high-spin *ferro* configurations.
        # The current implementation assumes that magnetic species are s^2 d^n p^0!
        # loops over high and low spin configurations.
        for spin in ["low", "high"]:
          spinjob = oxjob / (spin + "-spin")
          spinjob.jobparams["structure"] = structure
          spinjob.jobparams["nelect"]    = -oxidation
          spinjob.jobparams["magmom"]    = spin, neighbors
          spinjob.vasp = input.vasp

string = cPickle.dumps(jobs.current)
reloaded = cPickle.loads(string)
for job, name in reloaded.walk_through(outdir="results"):
  print name
  job.compute(norun=True, repat=files.input, outdir=name)
# jobs.current.update(jobs.load())
# for job in jobs.walk_through():
#   print job.job["outdir"]
#   if "nelect" in job.job: print job.job["nelect"]



