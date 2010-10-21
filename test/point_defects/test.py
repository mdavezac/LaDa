""" Submodule defining the *third wave* to compute point-defects. """
__docformat__ = "restructuredtext en"
__all__ = ['pointdefect_wave']

def pointdefect_wave(path=None, inputpath=None, **kwargs):
  """ Creates point-defect wave using ground-state job-dictionary. 

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

      Creates a point-defect wave from the materials computed in the
      magnetic and non-magnetic waves. Usage is fairly simple. If the pickle
      for the magnetic/non-magnetic wave is called ``magnetic_wave``, then one
      need only open it and call the ``pointdefect_wave``.

      >>> explore all magnetic_wave 
      >>> import test
      >>> test.magnetic_wave()
      >>> launch scattered

      The above will add point-defect calculations for all meterials and
      lattices of the ``magnetic_wave`` job-dictionary and save it (to the same 
      path unless an argument is provided to ``magnetic_wave``). Note that
      changing the location of the current job-dictionary has no effect. It
      would be possible but sounds too error prone:

      >>> explore all magnetic_wave 
      >>> goto /Fe2AlO4 # has no effect on what point-defects are added.
      >>> goto /Al2FeO4 # has no effect on what point-defects are added.
      >>> import test
      >>> # creates point-defects for both Fe2AlO4 and Al2FeO4: location does not matter.
      >>> test.magnetic_wave() 
      >>> launch scattered

      Point-defects calculations are added to
      a material if and only if all existing magnetic/non-magnetic jobs for
      that material have completed successfully. Furthermore, only untagged
      materials are accepted. Hence, to disable Fe2AlO4 lattices from having
      point-defects added to it, one can simply tag it:

      >>> explore all magnetic_wave 
      >>> current_jobdict["/Fe2AlO4"].tag()
      >>> import test
      >>> test.magnetic_wave()
      >>> launch scattered

      Similarly to make sure point-defect calculations are *not* created for
      the b5 lattice of the Fe2AlO4 material, one could tag it as follows:

      >>> explore all magnetic_wave 
      >>> current_jobdict["/Fe2AlO4/b5"].tag()
      >>> import test
      >>> test.magnetic_wave()
      >>> launch scattered
  """
  from tempfile import NamedTemporaryFile
  from os.path import dirname, normpath, relpath, join
  from IPython.ipapi import get as get_ipy
  from numpy import array, sum, abs
  from lada.jobs import JobDict
  from lada.vasp import read_input
  from lada.opt import Input
  from lada.crystal import point_defects as ptd

  # Loads jobdictionary and path as requested. 
  ip = get_ipy()
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
      input.update(read_input(file.name))
  if hasattr(jobdict, "maginput"):
    with NamedTemporaryFile() as file: 
      file.write(jobdict.nonmaginput)
      file.flush()
      input.update(read_input(file.name))
  if inputpath != None:
    input.update(read_input(inputpath))
    with open(inputpath, "r") as file: jobdict.maginput = file.read()
  input.update(kwargs)

  # saves inputfile to jobdictioanry if needed.
  if inputpath != None:
    input.update(read_input(inputpath))
    with open(inputpath, "r") as file: jobdict.pointdefectinput = file.read()
  # saves current script tof file.
  with open(__file__, "r") as file: jobdict.pointdefectscript = file.read()
  
  nb_new_jobs = 0
  # loops over completed structural jobs.
  for name in magnetic_groundstates():
    # gets the ground-states job-dictionary.
    groundstate = jobdict[name]
    # checks that the lattice and material are not tagged. 
    if groundstate[".."].is_tagged: continue
    if groundstate["../.."].is_tagged: continue
    # extracts the structure from it.
    superstructure, lattice = create_superstructure(groundstate, input)
    # extracts material.
    material = groundstate.material
    # extracts description of species.
    species = groundstate.functional.vasp.species

    # loop over substitutees.
    for B, substituters in input.point_defects.items():
      # loop over subtituters.
      for A in substituters:
        # loop over inequivalent point-defects sites.
        for structure, defect in ptd.all_defects(superstructure, lattice, B, A):
          # loop over oxidations states.
          for nb_extrae, oxname in ptd.charged_states(species, A, B):
            
            # creates list of moments. 
            new_moments = deduce_moment(A, species) 
            if len(new_moments) > 1: 
              moments = [ (min(new_moments), "/ls"), (max(new_moments), "/hs") ]
            else:
              moments = [ (max(new_moments), "") ]
            # loop  over moments.
            for moment, suffix in moments:
              name =  "PointDefects/{1}/{2}{0}".format(suffix, structure.name, oxname)
              
              # checks if job already exists. Does not change job if it exists!
              if name in groundstate[".."]: continue

              # creates new job.
              jobdict = groundstate["../"] / name
              jobdict.functional = input.relaxer
              jobdict.jobparams  = groundstate.jobparams.copy()
              jobdict.jobparams["structure"] = structure.deepcopy()
              jobdict.jobparams["nelect"] = nb_extrae
              jobdict.jobparams["relaxation"] = "ionic"
              jobdict.jobparams["ispin"] = 2
              jobdict.jobparams["set_symmetries"] = "off"
              jobdict.lattice  = lattice
              jobdict.material = material
              jobdict.defect   = defect
              # adds, modifies, or remove moment depending on defect type.
              if hasattr(superstructure, "magmom") or abs(moment) > 1e-12: 
                jstruct = jobdict.jobparams["structure"]
                # construct initial magmom
                if hasattr(superstructure, "magmom"):
                  jstruct.magmom = [u for u in superstructure.magmom]
                else: 
                  jstruct.magmom = [0 for u in superstructure.atoms]
                # now modifies according to structure.
                if B == None: # interstitial:
                  jstruct.magmom.append(moment)
                elif A == None: # vacancy -> remove moment.
                  jstruct.magmom.pop(defect.index)
                else: 
                  jstruct.magmom[defect.index] = moment
                # only keep moment if there are moments. 
                if sum(abs(jstruct.magmom)) < 1e-12 * float(len(jstruct.atoms)): del jstruct.magmom

              nb_new_jobs += 1

  # now saves new job dictionary
  print "Created {0} new jobs.".format(nb_new_jobs)
  if nb_new_jobs == 0: return
  ip.user_ns["current_jobdict"] = jobdict.root
  ip.magic("savejobs " + path)
            
          

def create_superstructure(groundstate, input):
  """ Creates a superstructure from existing structure. """
  from os.path import dirname, join
  from operator import itemgetter
  from numpy import dot
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
  extract = groundstate.functional.Extract( join(rootdir, groundstate.name[1:]) )
  assert extract.success, RuntimeError("Ground-state computation was not successful.")
  lattice = extract.structure.to_lattice()

  # creates superstructure.
  cell = dot(lattice.cell, input.supercell)
  result = fill_structure(cell, lattice)
  assert len(result.atoms) != len(lattice.sites), \
         ValueError("Superstructure as large as lattice. Disable this line if that's ok.")

  # adds magnetic moment if necessary.
  if hasattr(orig_lattice, "magmom"):
    assert extract.magnetization.shape[0] == len(lattice.sites),\
           RuntimeError("Could not find magnetization in ground-state's OUTCAR.")
    mlat = lattice.deepcopy()
    for atom, m in zip(mlat.sites, extract.magnetization[:,-1]):
      if abs(m) < 0.1: atom.type = '0'
      elif m < 0e0: atom.type = str(int(m-1))
      else: atom.type = str(int(m+1))
    moments = fill_structure(cell, mlat)
    result.magmom = [ int(i.type) for i in moments.atoms ]

  return result, lattice

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
  for nonmag in collect.grep(".*/.*/non-magnetic"):
    # check for success of all jobs (except for Point-defects).
    success = [u[1] for u in nonmag["../"].success.items() if u[0].find("PointDefects") == -1]
    if not all(success): continue
    # checks for lowest energy structure.
    energies = [u for u in nonmag["../"].total_energy.items() if u[0].find("PointDefects") == -1]
    energies = sorted(energies, key=itemgetter(1))
    yield energies[0][0]

def deduce_moment(type, species):
  """ Returns moment.

      This is a helper function which all atomic species the same with respect
      to the attribute ``moment``. If specie has no ``moment`` attribute,
      returns ``[0]``. If it exists and is a scalar, returns ``[moment]``. And
      if already is a list, returns as is.
  """
  if type == None: return [0]
  if not isinstance(type, str): type = type[0]
  if not hasattr(species[type], "moment"): return [0]
  if not hasattr(species[type].moment, "__iter__"):
    return [species[type].moment]
  return species[type].moment
