# will be params:
def create_jobdict(filename="input.py"):
  """ Returns a jobdictionary of point-defects.
     
      Parameters are set in input.py.
  """
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
  input = read_input(filename, input_dict)
  # sets as default lattice.
  input.lattice.set_as_crystal_lattice()

  # creates job dictionary.
  jobdict = jobs.JobDict()

  # creates super-structure
  supercell = fill_structure(input.supercell)

  # loop substitutions A on B. 
  for B, substituters in input.substitutions.items():
    for A in substituters:
      # loop over inequivalent substitution sites.
      for structure, substitution in ptd.substitution(supercell, input.lattice, B, A):

        # loop over oxidations states.
        for oxidation, name in ptd.charged_states(input.vasp.species, B, A):
          oxjob = jobdict / structure.name / name
 
          # now finds first neighbors. 12 is the highest coordination number, so
          # this should include the first shell.
          neighbors = [n for n in Neighbors(structure, 12, substitution.pos)]
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
              spinjob.jobparams["structure"] = structure
              spinjob.jobparams["nelect"]    = oxidation
              spinjob.jobparams["magmom"]    = spin, neighbors
              spinjob.functional = input.vasp

  return jobdict
