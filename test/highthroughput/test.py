#!/usr/bin/env python
""" High-Thoughput of A2BO4 structures. """
__docformat__ = "restructuredtext en"

def nonmagnetic_wave(path, inputpath="input.py", **kwargs):
  """ Jobs to explore possible ground-states. 
  
      :Parameters:
        path 
          Path where the job-dictionary will be saved. Calculations will be
          performed in the parent directory of this file. Calculations will be
          performed in same directory as this file.
        inputpath
          Path to an input file. Defaults to input.py. 
        kwargs
          Any keyword/value pair to add to the job-dictionaries' jobparams. These
          values will take precedence over anything in the input file.

      Creates a high-throughput job-dictionary to compute the non-magnetic
      ground-state of a host-material.  The new job-dictionary is loaded into
      memory automatically. No need to call explore. It is save to the path
      provided on input.
  """
  from re import compile
  import IPython.ipapi
  from lada.vasp import read_input
  from lada.jobs import JobDict
  from lada.crystal import fill_structure

  # reads input.
  input = read_input(inputpath)

  # sanity checks.
  for lattice in input.lattices:
    assert len(lattice.name) != 0, ValueError("Lattice has no name.")
  
  # regex
  specie_regex = compile("([A-Z][a-z]?)2([A-Z][a-z]?)([A-Z][a-z]?)4")

  # Job dictionary.
  jobdict = JobDict()
  if hasattr(input, "nbantiferro"): jobdict.nbantiferro = input.nbantiferro

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

    # loop over lattices. 
    for lattice in input.lattices:

      # creates a structure.
      structure = fill_structure(lattice.cell, lattice)
      # changes atomic species.
      for atom in structure.atoms:  atom.type  = species_dict[atom.type]
      # assigns it a name.
      structure.name = "{0} in {1}, spin-unpolarized.".format(material, lattice.name)
      # gets its scale.
      structure.scale = input.scale(structure)
  
      # job dictionary for this lattice.
      lat_jobdict = jobdict / material 
  
      job = lat_jobdict / lattice.name / "non-magnetic"
      job.functional = input.relaxer
      job.jobparams["structure"] = structure
      job.jobparams["ispin"] = 1
      # saves some stuff for future reference.
      job.material = material
      job.lattice  = lattice


  ip = IPython.ipapi.get()
  ip.user_ns["current_jobdict"] = jobdict
  ip.magic("savejobs " + path)


def magnetic_wave(path=None, jobdict=None, inputpath=None, nbantiferro=None, **kwargs):
  """ Creates magnetic wave from knowledge of previous wave. 

      :Parameters:
        path : str or None
          Path where the job-dictionary will be saved. Calculations will be
          performed in the parent directory of this file. If None, will use the
          current job-dictionary path.
        jobdict : `lada.jobs.JobDict` or None
          A job-dictionary containing the non-magnetic wave. If None, will
          refer to currently loaded job-dictionary.
        inputpath : str or None
          Path to an input file. If not present, then no input file is read and
          all parameters are taken from the non-magnetic wave.
        nbantiferro : int or None
          Number of random anti-ferro runs. If absent, looks into the input
          file given by inputpath. If that is absent as well, defaults to a variable located at
          the root of the jobdict (``jobdict.root.nbantiferro``), generally
          read from the input file at the time the non-magnetic wave was
          created. If that is absent as well (eg, wasn't in original
          non-magnetic input file), defaults to 0.
        kwargs
          Any keyword/value pair to add to the job-dictionaries' jobparams. These
          values will take precedence over anything in the input file.

      Creates magnetic wave from pre-existing non-magnetic wave. If no input
      file is given on input, then all parameters are obtained from the
      corresponding non-magnetic wave. 

      The new job-dictionary is loaded into memory automatically. No need to
      call explore. It is save to the path provided on input (or to the current
      job-dictionary path ifnot provided).  It will contain magnetic and
      non-magnetic calculations both. Pre-existing magnetic
      calculations will *not* be overwritten. However, additional anti-ferro
      configurations can be calculated by giving a large enough ``nbantiferro``
      on input.
  """
  from os.path import dirname, normpath, relpath, join
  import IPython.ipapi
  from lada.jobs import JobDict
  from lada.vasp import read_input

  ip = IPython.ipapi.get()
  if jobdict == None:
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
  else: basedir = dirname(path)
      
  # reads input.
  if inputpath != None:
    input = read_input(inputpath)
    if nbantiferro == None and hasattr(input, "antiferro"):
      nbantiferro = input.nbantiferro
  if nbantiferro == None: nbantiferro = 0

  nonmagname = "non-magnetic"
  has_changed = False
  for nonmagjob, name in jobdict.walk_through():
    # avoid other jobs (eg magnetic jobs).
    basename = normpath("/" + name + "/../")
    if relpath(name, basename[1:]) != nonmagname: continue
    # avoid tagged jobs.
    if nonmagjob.is_tagged: continue
    # check for success and avoid failures.
    extract = nonmagjob.functional.Extract(join(basedir, name)) 
    if not extract.success: continue
    if not is_magnetic_system(extract.structure, extract.functional.species): continue

    # now tries and creates high-spin ferro jobs if it does not already exist.
    jobname = normpath(basename + "/hs_ferro")
    magmom = hs_ferro(extract.structure, extract.functional.species)
    if magmom != None and jobname not in jobdict:
      job = jobdict / jobname
      job.functional = input.relaxer if inputpath != None else nonmagjob.functional
      job.jobparams["structure"] = extract.structure
      job.jobparams["structure"].name = "{0} in {1}, high-spin.".format(material, lattice.name)
      job.jobparams["magmom"] = magmom
      job.jobparams["ispin"] =  2
      job.jobparams.update(kwargs)
      # saves some stuff for future reference.
      job.material = nonmagjob.material
      job.lattice  = nonmagjob.lattice
      has_changed = True

    # now tries and creates low-spin ferro jobs if it does not already exist.
    jobname = normpath(basename + "/ls_ferro")
    magmom = ls_ferro(extract.structure, extract.functional.species)
    if magmom != None and jobname not in jobdict:
      job = jobdict / jobname
      job.functional = input.relaxer if inputpath != None else nonmagjob.functional
      job.jobparams["structure"] = extract.structure
      job.jobparams["structure"].name = "{0} in {1}, low-spin.".format(material, lattice.name)
      job.jobparams["magmom"] = magmom
      job.jobparams["ispin"] =  2
      job.jobparams.update(kwargs)
      # saves some stuff for future reference.
      job.material = nonmagjob.material
      job.lattice  = nonmagjob.lattice
      has_changed = True


    # now tries and creates anti-ferro-lattices jobs if it does not already exist.
    magmom = sublatt_antiferro(extract.structure, extract.functional.species) 
    jobname = normpath(basename + "/anti-ferro-0")
    if magmom != None and jobname not in jobdict:
      job = jobdict / jobname
      job.functional = input.relaxer if inputpath != None else nonmagjob.functional
      job.jobparams["structure"] = extract.structure
      job.jobparams["structure"].name = "{0} in {1}, lattice anti-ferro.".format(material, lattice.name)
      job.jobparams["magmom"] = magmom
      job.jobparams["ispin"] =  2
      job.jobparams.update(kwargs)
      # saves some stuff for future reference.
      job.material = nonmagjob.material
      job.lattice  = nonmagjob.lattice
      has_changed = True

    # random anti-ferro.
    for i in range(1, 1+nbantiferro):
      magmom = random(extract.structure, extract.functional.species)
      if magmom == None: continue
      jobname = normpath("/" + basename + "/anti-ferro-{0}".format(i))
      if jobname in jobdict: continue
      job = jobdict / jobname
      job.functional = input.relaxer if inputpath != None else nonmagjob.functional
      job.jobparams["structure"] = extract.structure
      job.jobparams["structure"].name = "{0} in {1}, random anti-ferro.".format(material, lattice.name)
      job.jobparams["magmom"] = magmom
      job.jobparams["ispin"] =  2
      job.jobparams.update(kwargs)
      # saves some stuff for future reference.
      job.material = nonmagjob.material
      job.lattice  = nonmagjob.lattice
      has_changed = True


  if not has_changed:
    print "No new jobs. "
    return
  ip = IPython.ipapi.get()
  ip.user_ns["current_jobdict"] = jobdict.root
  ip.magic("savejobs " + path)

def is_magnetic_system(structure, species):
  """ True if system is magnetic. """
  from lada.crystal import specie_list

  for u in [u for name, u in species.items() if name in specie_list(structure)]:
    if not hasattr(u, "moment"): continue
    if not hasattr(u.moment, "__iter__"): 
      if abs(u) > 1e-12: return True
      continue
    for a in u.moment:
      if abs(a) > 1e-12: return True
    
  return False

def ls_ferro(structure, species):
  """ Returns magmom VASP flag for low-spin ferromagnetic order. """
  from lada.crystal import specie_list

  magmom = ""
  for specie in specie_list(structure):
    assert specie in species,\
           KeyError("specie {0} not found in pseudo-potential dictionary.".format(specie))
    atoms = [atom for atom in structure.atoms if atom.type == specie]
    if hasattr(species[specie], "moment"):
      moment = species[specie].moment
      if hasattr(moment, "__iter__"): moment = min(moment)
      magmom += "{0}*{1:.2f} ".format(len(atoms), moment)
    else: magmom += "{0}*0 ".format(len(atoms), 0)
  return magmom

def hs_ferro(structure, species):
  """ Returns magmom VASP flag for high-spin ferromagnetic order. """
  from lada.crystal import specie_list

  magmom, has_both = "", False
  for specie in specie_list(structure):
    assert specie in species,\
           KeyError("specie {0} not found in pseudo-potential dictionary.".format(specie))
    atoms = [atom for atom in structure.atoms if atom.type == specie]
    if hasattr(species[specie], "moment"):
      moment = species[specie].moment
      if hasattr(moment, "__iter__"):
        moment = max(moment)
        has_both = True
      magmom += "{0}*{1:.2f} ".format(len(atoms), moment)
    else: magmom += "{0}*0 ".format(len(atoms), 0)
  return magmom if has_both else None

def sublatt_antiferro(structure, species):
  """ Anti ferro order with each cation type in a different direction. """
  from lada.crystal import specie_list

  magmom, sign, nb = "", 1e0, 0
  for specie in specie_list(structure):
    assert specie in species,\
           KeyError("specie {0} not found in pseudo-potential dictionary.".format(specie))
    atoms = [atom for atom in structure.atoms if atom.type == specie]
    if hasattr(species[specie], "moment"):
      moment = species[specie].moment
      if hasattr(moment, "__iter__"): moment = max(moment)
      magmom += "{0}*{1:.2f} ".format(len(atoms), moment * sign)
      nb += 1
      sign *= -1e0
    else: magmom += "{0}*0 ".format(len(atoms), 0)
  return None if nb < 2 else magmom

def random(structure, species):
  """ Random magnetic order. """
  from random import uniform
  from lada.crystal import specie_list

  magmom = ""
  for specie in specie_list(structure):
    assert specie in species,\
           KeyError("specie {0} not found in pseudo-potential dictionary.".format(specie))
    atoms = [atom for atom in structure.atoms if atom.type == specie]
    ps = species[specie]
    if hasattr(species[specie], "moment"): 
      moment = species[specie].moment
      if hasattr(moment, "__iter__"): moment = max(moment)
      for atom in atoms: magmom += "{0:.2f} ".format( uniform(-moment, moment) )
    else: magmom += "{0}*0   ".format(len(atoms))

  return magmom

