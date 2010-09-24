def ground_state_wave(path, inputpath="input.py"):
  """ Jobs to explore possible ground-states. 
  
      :Param path: 
        Path where the job-dictionary will be saved. Calculations will be
        performed in the parent directory of this file.
      :Param input: 
        Path to an input file. Defaults to input.py

      Creates a high-throughput job-dictionary to compute the magnetic
      ground-state of a host-material.
  """
  import IPython.ipapi
  from lada.vasp import read_input
  from lada.jobs import save

  # reads input.
  input = read_input(inputpath)

  # sanity checks.
  assert len(input.lattice.name) != 0, ValueError("Lattice has no name.")
  
  # regex
  specie_regex = compile("([A-Z][a-z]?)2([A-Z][a-z]?)([A-Z][a-z]?)4")

  # Job dictionary.
  jobdict = jobs.JobDict()

  # loop over materials.
  for material in input.materials:

    # creates dictionary to replace A2BX4 with meaningfull species.
    match = specie_regex.match(material)
    assert match != None, RuntimeError("Incorrect material " + material + ".")
    # checks species are known to vasp functional
    for i in range(1, 4):
      assert match.group(i) in input.vasp.species,\
             RuntimeError("%s not in specie dictionary." % (match.group(i)))
    # actually creates dictionary.
    species_dict = {"A": match.group(1), "B": match.group(2), "X": match.group(3)}

    # creates a structure.
    structure = fill_structure(input.lattice.cell, input.lattice)
    # assigns it a name.
    structure.name = "{0} in {1}, spin-unpolarized.".format(material, input.lattice.name, 
    # gets its scale.
    structure.scale = input.scale(structure)
    # changes atomic species.
    for atom in structure.atoms:  atom.type  = species_dict[atom.type]

    # job dictionary for this lattice.
    lat_jobdict = jobdict / material 

    job = lat_jobdict / "non-magnetic"
    job.functional = input.relaxer
    job.jobparams["structure"] = structure
    job.jobparams["ispin"] = 1

  ip = IPython.ipapi.get()
  ip.user_ns["current_jobdict"] = jobdict
  ip.magic("savejob " + path)



# will be params:
def _ox_spin_loops(A, B, defect, structure, input, is_interstitial=False):
  """ Loops over oxidation and moments. """
  from copy import deepcopy
  from lada.crystal import fill_structure, Neighbors, point_defects as ptd
  from lada.jobs import JobDict
  from lada.vasp.methods import RelaxIons

  # creates jobdictionary 
  jobdict = JobDict()
  oxidation = input.vasp.species[defect.type].oxidation if defect.type != "None" else 0
  # loop over oxidations states.
  for nb_extrae, oxname in ptd.charged_states(input.vasp.species, A, B):

    if is_interstitial: nb_extrae *= -1
    moment = nb_extrae - oxidation
    # loop over low spin magnetic states, integer and average.
    iterspins = ptd.low_spin_states(structure, defect, input.vasp.species, moment)
    for indices, moments in iterspins:

      job = jobdict / oxname / ptd.magname(moments, "moment")
      job.jobparams["structure"] = structure
      job.jobparams["nelect"] = nb_extrae
      job.jobparams["nupdown"] = sum(moments)
      job.jobparams["magmom"] = ptd.magmom(indices, moments, len(structure.atoms))
      job.functional = RelaxIons(input.vasp, first_trial = input.relaxation_parameters)
      job.defect = deepcopy(defect) # original atom prior to defect.

    # loop over high spin magnetic states, integer and average.
    iterspins = ptd.high_spin_states( structure, defect, input.vasp.species, moment)
    for indices, moments in iterspins:

      job = jobdict / oxname / ptd.magname(moments, "moment")
      job.jobparams["structure"] = structure
      job.jobparams["nelect"] = moment
      job.jobparams["nupdown"] = sum(moments)
      job.jobparams["magmom"] = ptd.magmom(indices, moments, len(structure.atoms))
      job.functional = RelaxIons(input.vasp, first_trial = input.relaxation_parameters)
      job.defect = deepcopy(defect) # original atom prior to defect.

    # do paramagnetic calculation.
    job = jobdict / oxname / "paramagnetic"
    job.jobparams["structure"] = structure
    job.jobparams["nelect"] = nb_extrae
    job.jobparams["nupdown"] = None
    job.jobparams["magmom"] = None
    job.functional = RelaxIons(input.vasp, first_trial = input.relaxation_parameters)
    job.defect = deepcopy(defect) # original atom prior to defect.

  return jobdict


def create_jobdict(filename="input.py"):
  """ Returns a jobdictionary of point-defects.
     
      Parameters are set in input.py.
  """
  import cPickle
  from sys import exit
  from os.path import join

  from lada.opt import read_input
  from lada.vasp import Vasp, Specie, specie, files
  from lada.crystal import fill_structure, Neighbors, point_defects as ptd
  from copy import deepcopy
  from lada import jobs

  # names we need to create input.
  input_dict = { "Vasp": Vasp, "U": specie.U, "nlep": specie.nlep }
  # reads input.
  input = read_input(filename, input_dict)
  # sets as default lattice.
  input.lattice.set_as_crystal_lattice()
  if not hasattr(input, "relaxation_parameters"): input.relaxation_parameters = {}

  # creates job dictionary.
  jobdict = jobs.JobDict()

  # creates super-structure
  supercell = fill_structure(input.supercell)

  # loop substitutions A on B. 
  for B, substituters in input.substitutions.items():
    for A in substituters:

      # loop over inequivalent substitution sites.
      for structure, substitution in ptd.substitution(supercell, input.lattice, B, A):

        # loops over oxidation and moments are collected in _ox_spin_loop,
        # since it can be used for intersitials.
        jobdict[structure.name] = _ox_spin_loops(A, B, substitution, structure, input)
  
  # return if not interstitials.
  if not hasattr(input, "interstitials"): return jobdict
  
  # loop over interstitials.
  for type, positions in input.interstitials.items():

    # loop over substitutional positions (and name).
    for position in positions: 

      # create structures.
      structure = fill_structure(input.supercell)
      structure.add_atom = position[:-1], type
      structure.name = "{0}_interstitial_{1}".format(type, position[3])

      # loops over oxidation and moments are collected in _ox_spin_loop,
      # since it can be used for substitutionals.
      defect = deepcopy(structure.atoms[-1])
      defect.type = "None"
      jobdict[structure.name] = _ox_spin_loops(type, None, defect, structure, input, True)
       
      

  return jobdict
