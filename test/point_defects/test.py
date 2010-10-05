def pointdefect_wave(path=None, inputpath=None, **kwargs):
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
  if path == None:
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

  nb_new_jobs = 0
  # loops over completed structural jobs.
  for name in magnetic_groundstates():
    # gets the ground-states job-dictionary.
    groundstate = jobdict[name]
    # extracts the structure from it.
    superstructure = create_superstructure(groundstate, input)

    # loop over substitution.
    for B, substituters in substitutions.items():
      for A in substituters:
        # loop over inequivalent substitution sites.
        for structure, substitution in ptd.substitution(superstructure, input.lattice, B, A):
          # gets moments.
          new_moment = deduce_moment(A, input.vasp.species)

          # low-spin or no-spin case.
          # creates new job-dictionary name.
          if len(new_moments) == 2:   name =  "PointDefects/ls-{0}".format(structure.name)
          elif len(new_moments) == 1: name = "PointDefects/{0}".format(structure.name)
          else: raise RuntimeError("Too many moments in specie description.")

          # checks if job already exists.
          if name not in groundstate[".."]: 
            jobdict = groundstate["../"] / name
            jobdict.functional = groundstates.functional
            jobdict.jobparams  = groundstates.jobparams.copy()
            jobdict.jobparams["structure"] = deepcopy(structure)
            jobdict.lattice  = groundstate.lattice
            jobdict.material = groundstate.material
            jobdict.defect   = substitution
            # adds or remove moment.
            if new_moments[0] == None: # vacancy -> remove moment.
              jobdict.jobparams["structure"].magmom.pop(substitution.index)
            else: 
              jobdict.jobparams["structure"].magmom[substitution.index] = min(new_moments)
            nb_new_jobs += 1

          # high-spin case.
          name =  "PointDefects/hs-{0}".format(structure.name)
          if len(new_moments) == 1 and name not in groundstate[".."]: 
            jobdict = groundstate["../"] / name
            jobdict.functional = groundstates.functional
            jobdict.jobparams  = groundstates.jobparams.copy()
            jobdict.jobparams["structure"] = deepcopy(structure)
            jobdict.lattice  = groundstate.lattice
            jobdict.material = groundstate.material
            jobdict.defect   = substitution
            # adds or remove moment.
            if new_moments[0] == None: # vacancy -> remove moment.
              jobdict.jobparams["structure"].magmom.pop(substitution.index)
            else: 
              jobdict.jobparams["structure"].magmom[substitution.index] = max(new_moments)
            nb_new_jobs += 1
          
    # loop over interstitials.
    for type, positions in input.interstitials.items():
      # loop over substitutional positions (and name).
      for position in positions: 
        # create structures.
        structure = superstructure.copy()
        structure.add_atom = position[:-1], type
        structure.name = "{0}_interstitial_{1}".format(type, position[3])
  
          # gets moments.
          new_moment = deduce_moment(A, input.vasp.species)

          # low-spin or no-spin case.
          # creates new job-dictionary name.
          if len(new_moments) == 2:   name =  "PointDefects/ls-{0}".format(structure.name)
          elif len(new_moments) == 1: name = "PointDefects/{0}".format(structure.name)
          else: raise RuntimeError("Too many moments in specie description.")

          # checks if job already exists.
          if name not in groundstate[".."]: 
            jobdict = groundstate["../"] / name
            jobdict.functional = groundstates.functional
            jobdict.jobparams  = groundstates.jobparams.copy()
            jobdict.jobparams["structure"] = deepcopy(structure)
            jobdict.lattice  = groundstate.lattice
            jobdict.material = groundstate.material
            jobdict.defect   = substitution
            # adds or remove moment.
            if new_moments[0] == None: # vacancy -> remove moment.
              jobdict.jobparams["structure"].magmom.pop(substitution.index)
            else: 
              jobdict.jobparams["structure"].magmom[substitution.index] = min(new_moments)
            nb_new_jobs += 1

          # high-spin case.
          name =  "PointDefects/hs-{0}".format(structure.name)
          if len(new_moments) == 1 and name not in groundstate[".."]: 
            jobdict = groundstate["../"] / name
            jobdict.functional = groundstates.functional
            jobdict.jobparams  = groundstates.jobparams.copy()
            jobdict.jobparams["structure"] = deepcopy(structure)
            jobdict.lattice  = groundstate.lattice
            jobdict.material = groundstate.material
            jobdict.defect   = substitution
            # adds or remove moment.
            if new_moments[0] == None: # vacancy -> remove moment.
              jobdict.jobparams["structure"].magmom.pop(substitution.index)
            else: 
              jobdict.jobparams["structure"].magmom[substitution.index] = max(new_moments)
            nb_new_jobs += 1
        # loops over oxidation and moments are collected in _ox_spin_loop,
        # since it can be used for substitutionals.
        defect = deepcopy(structure.atoms[-1])
        defect.type = "None"
        dummy = (root / "PointDefects" / structure.name)
        has_changed |= dummy.add_new( charge_and_spins(A, B, defect, structure, input) ) 

def create_superstructure(groundstate, input):
  """ Creates a superstructure from existing structure. """
  from os.path import dirname
  from operator import itemgetter
  from numpy import matrix
  from IPython.ipapi import get as get_ipy
  from lada.crystal import fill_structure

  # sanity checks,
  assert "structure" in groundstate.jobparams,\
         ValueError("Could not find structure in ground-state job-dictionary.")
  assert hasattr(groundstate.functional, "Extract"),\
         ValueError("Could not find extraction class in ground-state job-dictionary.")
  
  ip = get_ipy()
  assert "current_jobdict_path" in ip.user_ns,\
         RuntimeError("Could not find path for current job-dictionary.")
  rootdir = dirname(ip.user_ns["current_jobdict_path"])
  # gets original lattice from job-dictionary.
  orig_lattice = groundstate.jobparams["structure"].to_lattice()
  
  # Extracts computed lattice from ground state calculation.
  extract = groundstate.functional.Extract( join(rootdir, groundstate.name) )
  assert extract.success, RuntimeError("Ground-state computation was not successful.")
  lattice = extract.structure.to_lattice()

  # creates superstructure.
  cell = matrix(lattice.cell, dtype="float64") * matrix(input.supercell, dtype="float64")
  result = fill_structure(cell, lattice)

  # adds magnetic moment if necessary.
  if hasattr(orig_lattice, "magmom"):
    assert extract.magnetization.shape[0] == len(lattice.cell),\
           RuntimeError("Could not find magnetization in ground-state's OUTCAR.")
    mlat = lattice.copy()
    for atom, m in zip(lattice.atoms, extract.magnetization[:,-1]):
      i = sorted( enumerate([abs(m), abs(m-1e0), abs(m+1e0), abs(m-5e0), abs(m+5e0)]),\
                  key=itemgetter(1))[0]
      atom.type = str(i)
    moments = fill_structure(cell, mlat)
    result.magmom = [ [0, 1, -1, 5, -5][int(i.type)] for i in moments.atoms ]

  return result

def magnetic_groundstates():
  """ Yields name of magnetic-groundstates from current job-dictionary.

      A set of magnetic-states for the same lattice and materials is defined by
      all jobs residing in the parent directory of .*/.*/non-magnetic, other
      than PointDefects. 

      All jobs within a set of magnetic-states must be
      finished. Otherwise, that particular combination of material + lattice is
      not considered. 

      This yields the fully qualified job-name of each lowest energy magnetic
      ground-state within the current job-dictionary.
  """
  from operator import itemgetter
  from lada.ipython import Collect
  collect = Collect()
  # loops over untagged non-magnetic structural jobs.
  for nonmag in collect.grep("/.*/.*/non-magnetic"):
    # check for success of all jobs (except for Point-defects).
    success = [u[1] for u in nonmag["../"].success.items() if u[0].find("PointDefects") == -1]
    if not all(success): continue
    # checks for lowest energy structure.
    energies = [u for u in nonmag["../"].total_energy.items() if u[0].find("PointDefects") == -1]
    energies = sorted(energies, key=itemgetter(1))
    yield energies[-1][0]

def deduce_moment(type, species):
  """ Returns moment.

      This is a helper function which all atomic species the same with respect
      to the attribute ``moment``. If specie has no ``moment`` attribute,
      returns ``[0]``. If it exists and is a scalar, returns ``[moment]``. And
      if already is a list, returns as is.
  """
  if type == None: return [None]
  if not hasattr(species[type], "moment"): return [0]
  if not hasattr(species[type].moment, "__iter__"):
    return [species[type].moment]
  return species[type].moment
