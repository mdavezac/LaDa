#!/usr/bin/env python
""" High-Thoughput of A2BO4 structures. """
__docformat__ = "restructuredtext en"
__all__ = ['nonmagnetic_wave', 'magnetic_wave']

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
          Any keyword/value pair to take precedence over anything in the input
          file.

      Creates a high-throughput job-dictionary to compute the non-magnetic
      ground-state of a host-material.  The new job-dictionary is loaded into
      memory automatically. No need to call explore. It is save to the path
      provided on input.
  """
  from re import compile
  from IPython.core.interactiveshell import InteractiveShell
  from copy import deepcopy
  from lada.vasp import read_input
  from lada.jobs import JobDict
  from lada.crystal import fill_structure
  from lada import interactive

  # reads input.
  input = read_input(inputpath)
  input.update(kwargs)

  # sanity checks.
  for lattice in input.lattices:
    assert len(getattr(lattice, 'name', '') != 0, ValueError("Lattice has no name.")
  
  # regex
  specie_regex = compile("([A-Z][a-z]?)([A-Z][a-z]?)([A-Z][a-z]?)")

  # Job dictionary.
  jobdict = JobDict()

  # loop over materials.
  for material in input.materials:

    # creates dictionary to replace A2BX4 with meaningfull species.
    match = specie_regex.match(material)
    print match.group(1)
    print match.group(2)
    print match.group(3)
    assert match is not None, RuntimeError("Incorrect material " + material + ".")
    # checks species are known to vasp functional
    for i in range(1, 4):
      assert match.group(i) in input.vasp.species,\
             RuntimeError("No pseudo-potential defined for {0}.".format(match.group(i)))
    # actually creates dictionary.
    species_dict = {"A": match.group(1), "B": match.group(2), "X": match.group(3)}

    # loop over lattices. 
    for lattice in input.lattices:

      # creates a structure.
      structure = deepcopy(lattice)
      # changes atomic species.
      for atom in structure.atoms:  atom.type  = species_dict[atom.type]
      structure = vasp_ordered(structure)
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


  interactive.jobdict = jobdict
  InteractiveShell.instance().magic("savejobs " + path)


def magnetic_wave(path=None, inputpath='input.py', **kwargs):
  """ Creates magnetic wave for current job-dictionary.

      :Parameters:
        path : str or None
          Path where the modified job-dictionary will be saved. Calculations will be
          performed in the parent directory of this file. If None, will use the
          current job-dictionary path.
        inputpath : str or None
          Path to an input file. If not present, then no input file is read and
          all parameters are taken from the non-magnetic wave.
        kwargs
          Any keyword/value pair to take precedence over anything in the input file.

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
  from tempfile import NamedTemporaryFile
  from os.path import dirname, normpath, relpath, join
  from IPython.core.interactiveshell import InteractiveShell
  from lada.jobs import JobDict
  from lada.vasp import read_input
  from lada.opt import Input

  # Loads jobdictionary and path as requested. 
  if interactive.jobdict is None: 
    print "No current job-dictionary."
    return
  if interactive.jobdict_path is None: 
    print "No path for current job-dictionary."
    return
  basedir = dirname(interactive.jobdict_path)
      
  input = read_input(inputpath)

  # will loop over all jobs, looking for *successfull* *non-magnetic* calculations. 
  # Only magnetic jobs which do NOT exist are added at that point.
  nonmagname = "non-magnetic"
  nb_new_jobs = 0
  for name, nonmagjob in jobdict.iteritems():
    # avoid tagged jobs.
    if nonmagjob.is_tagged: continue
    # avoid other jobs (eg magnetic jobs).
    basename = normpath("/" + name + "/../")
    if relpath(name, basename[1:]) != nonmagname: continue
    # check for success and avoid failures.
    extract = nonmagjob.functional.Extract(join(basedir, name)) 
    if not extract.success: continue
    if not is_magnetic_system(extract.structure, extract.functional.species): continue

    # loads lattice and material from non-magnetic job.
    material = nonmagjob.material
    lattice = nonmagjob.lattice

    # figures out whether we have both high and low spins. 
    if has_high_and_low(extract.structure, extract.functional.species):
          hnl = [(min, "ls-"), (max, "hs-")]
    else: hnl = [(min, "")] 
    # now loops over moments.
    for func, prefix in hnl: 
      # now tries and creates high-spin ferro jobs if it does not already exist.
      jobname = normpath("{0}/{1}ferro".format(basename, prefix))
      structure, magmom = ferro(extract.structure, extract.functional.species, func)
      if magmom and jobname not in jobdict:
        structure.name = "{0} in {1}, {2}ferro."\
                         .format(material, lattice.name, prefix)
        job = jobdict / jobname
        job.functional = input.relaxer 
        job.jobparams["structure"] = structure.copy()
        job.jobparams["magmom"] = True
        job.jobparams["ispin"] =  2
        # saves some stuff for future reference.
        job.material = material
        job.lattice  = lattice
        nb_new_jobs += 1

      # now tries and creates anti-ferro-lattices jobs if it does not already exist.
      structure, magmom = species_antiferro(extract.structure, extract.functional.species, func) 
      jobname = normpath("{0}/{1}anti-ferro-0".format(basename, prefix))
      if magmom and jobname not in jobdict:
        structure.name = "{0} in {1}, {2}specie-anti-ferro."\
                         .format(material, lattice.name, prefix)

        job = jobdict / jobname
        job.functional = input.relaxer
        job.jobparams["structure"] = structure.copy()
        job.jobparams["magmom"] = True
        job.jobparams["ispin"] =  2
        # saves some stuff for future reference.
        job.material = material
        job.lattice  = lattice
        nb_new_jobs += 1

      # random anti-ferro.
      for i in range(1, 1+input.nbantiferro):
        structure, magmom = random(extract.structure, extract.functional.species, func)
        if not magmom: continue
        jobname = normpath("{0}/{1}anti-ferro-{2}".format(basename, prefix, i))
        if jobname in jobdict: continue
        structure.name = "{0} in {1}, random anti-ferro.".format(material, lattice.name)

        job = jobdict / jobname
        job.functional = input.relaxer if inputpath is not None else nonmagjob.functional
        job.jobparams["structure"] = structure.copy()
        job.jobparams["magmom"] = True
        job.jobparams["ispin"] =  2
        # saves some stuff for future reference.
        job.material = material
        job.lattice  = lattice
        nb_new_jobs += 1

  # now saves new job dictionary
  print "Created {0} new jobs.".format(nb_new_jobs)
  if nb_new_jobs == 0: return
  interactive.jobdict = jobdict.root
  InteractiveShell.instance().magic("savejobs " + path)

def is_magnetic_system(structure, species):
  """ True if system is magnetic. """
  from lada.crystal import specie_list

  for u in [u for name, u in species.items() if name in specie_list(structure)]:
    if not hasattr(u, "moment"): continue
    if not hasattr(u.moment, "__iter__"): 
      if abs(u.moment) > 1e-12: return True
      continue
    for a in u.moment:
      if abs(a) > 1e-12: return True
    
  return False

def has_high_and_low(structure, species):
  """ True if some species have both high and low spins. """
  for atom in structure.atoms:
    if len(deduce_moment(atom, species)) > 1: return True
  return False

def deduce_moment(atom, species):
  """ Returns moment.

      This is a helper function which all atomic species the same with respect
      to the attribute ``moment``. If specie has no ``moment`` attribute,
      returns ``[0]``. If it exists and is a scalar, returns ``[moment]``. And
      if already is a list, returns as is.
  """
  if not hasattr(species[atom.type], "moment"): return [0]
  if not hasattr(species[atom.type].moment, "__iter__"):
    return [species[atom.type].moment]
  return species[atom.type].moment

def ferro(structure, species, func=min):
  """ Returns magmom VASP flag for low-spin ferromagnetic order. """
  return [ func(deduce_moment(a, species)) for a in structure.atoms ]

def species_antiferro(structure, species, func=min):
  """ Low spin anti ferro order with each cation specie in a different direction. """
  from numpy import math, abs
  # checks whether lattice-sites are magnetic or not.
  result = structure.copy()
  signs = {}
  for atom in result:
    m = func(deduce_moment(atom, species))
    signs[atom.type] = 1 if abs(m) < 1e-12 else -1
  if len([u[1] for u in signs.items() if u[1] == -1]) < 2: return None
  # makes alternating sign list.
  dummy = 1
  for k in sorted(signs.keys()):
    if signs[k] != -1: continue
    signs[k] = 1 if dummy == 1 else -1
    dummy = -1 * dummy

  # adds magmom tag to structure's atoms.
  s = 0e0
  for atom in result:
    atom.magmom = float(signs[atom.type]) * func(deduce_moment(atom, species))
    s += abs(atom.magmom)
  return s > 1e-6, result

def random(structure, species, func=min):
  """ High-spin random magnetic order. """
  from random import choice
  return [ choice([-1e0, 1e0]) * func(deduce_moment(a, species))\
           for a in structure.atoms ]
