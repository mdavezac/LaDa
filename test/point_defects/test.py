def pointdefect_wave(inputpath=None, **kwargs):
  """ Creates point-defect wave using ground-state job-dictionary. """
  from tempfile import NamedTemporaryFile
  from os.path import dirname, normpath, relpath, join
  import IPython.ipapi
  from lada.jobs import JobDict
  from lada.vasp import read_input
  from lada.opt import Input

  # Loads jobdictionary and path as requested. 
  ip = IPython.ipapi.get()
  if "current_jobdict" not in ip.user_ns: 
    print "No current job-dictionary." 
    return
  jobdict = ip.user_ns["current_jobdict"].root
  if "current_jobdict_path" not in ip.user_ns:
    print "No known path for current dictionary and no path specified on input."
    return
  path = ip.user_ns["current_jobdict_path"]
  basedir = dirname(path)

  # create input dictionary. First reads non-magnetic input, then magnetic
  # input, then kwargs. Order of precedence is exact opposite.
  input = Input()
  if hasattr(jobdict, "nonmaginput"):
    with NamedTemporaryFile() as file: 
      file.write(jobdict.nonmaginput)
      file.flush()
      input.update(read_input(file.name).__dict__)
  if hasattr(jobdict, "maginput"):
    with NamedTemporaryFile() as file: 
      file.write(jobdict.nonmaginput)
      file.flush()
      input.update(read_input(file.name).__dict__)
  if inputpath != None:
    input.update(read_input(inputpath))
    with open(inputpath, "r") as file: jobdict.maginput = file.read()
  input.update(kwargs)
  
  assert hasattr(input, "supercell"), RuntimeError("Supercell not given on input.")
  assert hasattr(input, "interstitials") or hasattr(input, "substitutions"),\
         RuntimeError("Neither substitutions nor interstitials given on input.")
  interstitials = {} if not hasattr(input, "interstitials") else input.interstitials
  substitutions = {} if not hasattr(input, "substitutions") else input.substitutions

  # loops over completed structural jobs.
  for name in completed_structurals():
    root = jobdict[name]
    if not hasattr(root["non-magnetic"], "lattice"): continue
    # starts creating Point-Defect jobs.
    supercell = fill_structure(input.supercell, root["non-magnetic"].lattice)

    # loop over interstitials.
    for B, substituters in substitutions.items():
      for A in substituters:
        # loop over inequivalent substitution sites.
        for structure, substitution in ptd.substitution(supercell, input.lattice, B, A):
          # Add new jobs only.
          dummy = (root / "PointDefects" / structure.name)
          has_changed |= dummy.add_new( charge_and_spins(A, B, substitution, structure, input) ) 

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
        dummy = (root / "PointDefects" / structure.name)
        has_changed |= dummy.add_new( charge_and_spins(A, B, substitution, structure, input) ) 



  

def completed_structurals():
  """ Yields structural jobs which are complete.

      Returns the name of the job containing all magnetic jobs for each lattice
      and material if and only if all magnetic jobs are finished and
      successful. 
  """
  from lada.ipython import Collect
  collect = Collect()
  # loops over untagged non-magnetic structural jobs.
  for nonmag in collect.grep("/.*/.*/non-magnetic"):
    successes = collect["../"].success.items()
    if all( [value for key, value in successes if value.find("PointDefects") == -1] ):
      yield collect["../"].position








def charges_and_spins(A, B, defect, structure, input, is_interstitial=False):
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


def point_defects(jobdict, input)
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
