# will be params:
def create_jobdict():
  import cPickle
  from sys import exit
  from os.path import join

  from lada.opt import read_input
  from lada.vasp import Vasp, Specie, specie, files
  from lada.vasp.methods import RelaxCellShape
  from lada.crystal import fill_structure, Neighbors, point_defects as ptd
  from lada import jobs

  # names we need to create input.
  input_dict = { "Vasp": Vasp, "U": specie.U, "nlep": specie.nlep }
  # reads input.
  input = read_input("input.py", input_dict)
  # sets as default lattice.
  input.lattice.set_as_crystal_lattice()

  # creates job dictionary.
  jobdict = jobs.JobDict()

  # adds pure host calculation.
  jobdict.functional = RelaxCellShape(input.vasp, relaxation="volume cellshape ions")
  jobdict.jobparams["structure"] = fill_structure(input.lattice.cell)
  jobdict.jobparams["kpoints"] = "Automatic generation\n0\nMonkhorst\n8 8 8\n0 0 0"
  jobdict.jobparams["encut"] = 1.3
  jobdict.jobparams["set_symmetries"] = 1e-5


  # creates super-structure
  supercell = fill_structure(input.supercell)

  # loop over specie vacancy
  for A, specie in input.vasp.species.items(): 
    vacjobs = jobdict / ("vacancy_" + A)
    # loop over types of vacancies for given specie (eg symmetrically inequivalent sites).
    for structure, vacancy, name in ptd.vacancy(supercell, input.lattice, A):
      vacjob = vacjobs / name

      # loop over oxidation states.
      for charge, name in ptd.charged_states(specie):
        oxjob = vacjob / name

        # now finds first neighbors. 12 is the highest coordination number, so
        # this should include the first shell.
        neighbors = [n for n in Neighbors(structure, 12, vacancy.pos)]
        # only take the first shell and keep indices (to atom in structure) only.
        neighbors = [n.index for n in neighbors if n.distance < neighbors[0].distance + 1e-1]
        # reduce to those which are magnetic.
        neighbors = [n for n in neighbors if input.vasp.species[ structure.atoms[n].type ].magnetic]

        if len(neighbors) == 0: # no magnetic neighbors.
          oxjob.jobparams["structure"] = structure
          oxjob.jobparams["nelect"]    = charge
          oxjob.functional = input.vasp
          # Note: structure and nelect can now be accessed as
          # oxjob.structure and oxjob.nelect
        else: # has magnetic neighbors
          # Creates low spin and high-spin *ferro* configurations.
          # The current implementation assumes that magnetic species are s^2 d^n p^0!
          # loops over high and low spin configurations.
          for spin in ["low", "high"]:
            spinjob = oxjob / (spin + "-spin")
            spinjob.jobparams["structure"] =  structure
            spinjob.jobparams["nelect"]    =  charge
            spinjob.jobparams["magmom"]    =  (spin, neighbors)
            spinjob.functional = input.vasp
            # Note: structure, magmom, and nelect can now be accessed as
            # spinjob.structure, spinjob.magmom, and spinjob.nelect

  # loop substitutions.
  for A, B in [("Rh", "Zn"), ("Zn", "Rh") ]:
    # loop over inequivalent substitution sites.
    subjobs = jobdict / (A + "_on_" + B)
    nbineq = len(ptd.inequivalent_sites(input.lattice, A)) 
    for structure, substitution, name in ptd.substitution(supercell, input.lattice, A, B):
      subjob = subjobs / name

      # loop over oxidations.
      for oxidation, name in ptd.charged_states(input.vasp.species[A], input.vasp.species[B]):
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
          oxjob.jobparams["nelect"]    = oxidation
          oxjob.functional = input.vasp
        else: # has magnetic neighbors
          # Creates low spin and high-spin *ferro* configurations.
          # The current implementation assumes that magnetic species are s^2 d^n p^0!
          # loops over high and low spin configurations.
          for spin in ["low", "high"]:
            spinjob = oxjob / (spin + "-spin")
            spinjob["structure"] = structure
            spinjob["nelect"]    = oxidation
            spinjob["magmom"]    = spin, neighbors
            spinjob.functional = input.vasp

  return jobdict
