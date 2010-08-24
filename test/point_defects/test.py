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

        oxidation = input.vasp.species[substitution.type].oxidation 
        # loop over oxidations states.
        for nb_extrae, oxname in ptd.charged_states(input.vasp.species, B, A):

          moment = nb_extrae - oxidation
          # loop over low spin magnetic states, integer and average.
          iterspins = ptd.low_spin_states(structure, substitution, input.vasp.species, moment)
          for indices, moments in iterspins:

            job = jobdict / structure.name / oxname / ptd.magname(moments, "moment")
            job.jobparams["structure"] = structure
            job.jobparams["nelect"] = nb_extrae
            job.jobparams["nupdown"] = sum(moments)
            job.jobparams["magmom"] = ptd.magmom(indices, moments, len(structure.atoms))
            job.functional = input.vasp

          # loop over high spin magnetic states, integer and average.
          iterspins = ptd.high_spin_states( structure, substitution, input.vasp.species, moment)
          for indices, moments in iterspins:

            job = jobdict / structure.name / oxname / ptd.magname(moments, "moment")
            job.jobparams["structure"] = structure
            job.jobparams["nelect"] = moment
            job.jobparams["nupdown"] = sum(moments)
            job.jobparams["magmom"] = ptd.magmom(indices, moments, len(structure.atoms))
            job.functional = input.vasp

          # do paramagnetic calculation.
          job = jobdict / structure.name / oxname / "paramagnetic"
          job.jobparams["structure"] = structure
          job.jobparams["nelect"] = nb_extrae
          job.jobparams["nupdown"] = None
          job.jobparams["magmom"] = None
          job.functional = input.vasp


  return jobdict
