""" IPython interface for lada.
    ===========================

    IPython is the *enhanced interactive python shell*. In  practice this
    means it is a bash shell which incorporates interactive python programming,
    or the other way around. Think of it as a bash script where you don't need
    to call awk to add a few numbers, or print things in a different format,
    but rather use python's power to do everything natively and on the fly.
    It's a bash shell for the real geeks. It does have a `user guide`__, as
    well as a `tutorial`__ with plenty of tips.

    __ http://ipython.scipy.org/doc/stable/html/
    __ http://ipython.scipy.org/doc/manual/html/interactive/tutorial.html


    Installing the lada interface to ipython. 
    -----------------------------------------

    In order to initialize an ipython session with the lada interface, the
    following lines can be introduces in the ``main`` function of your
    ipy_user_conf.py

    >>> try: import lada.ipython 
    >>> except ImportError:
    >>>   print "Could not find module lada.ipython."
    >>>   pass
    >>> else: lada.ipython.ipy_init()

    This will give you access to all the lada magic functions. 

    In addition, I would recommend uncommenting/adding the following lines in
    you ipy_user_conf.py. They can be found at the beginning of the ``main``
    function.

    >>> import ipy_profile_sh

    This will allow you to use standard bash commands (cd, cp, find, cat, ...)
    without any fuss (Whithout this, bash commands can only be accessed by
    escaping them first with an exclamation mark, eg ``!cd ..`` or ``!rm -rf ~/\*``).

    >>> import ipy_editors
    >>> # Then for VI users.
    >>> ipy_editors.install_editor("vim +$line $file")    
    >>> # or for EMACS users.
    >>> ipy_editors.install_editor("emacs -nw +$line $file")

    There are many more things you can do. Checkout the two links given above.


    What is a Job-dictionary?
    -------------------------
    
    In practice a job-dictionary should be viewed as a directory tree where
    calculations will be performed. The old way was to create a bunch of
    directories with POSCAR, INCAR, POTCAR files at the start of the job.
    For example, if trying to compute Spinel and Olivine structures for the two
    materials ZnMgO and ZnRhO, the folling directory structure could have been created
    (by hand, or by stringing a bunch of bash scripts):

    :: 

      path/to/calcs/
        ZnRhO/
          Spinel/
            POSCAR  
            POTCAR
            INCAR 
          Olivine/
            POSCAR  
            POTCAR
            INCAR 
        ZnMgO/
          Spinel/
            POSCAR  
            POTCAR
            INCAR 
          Olivine/
            POSCAR  
            POTCAR
            INCAR 

    The game is then to launch calculations in each directory. The
    job-dictionary is a python object which replicates this type of tree-like
    structure. It's main advantage is that it can be manipulated
    programmatically within python,  with more ease than bash can manipulate a
    collections of files and directories. The programmatic approach (eg
    writing a script to create a job-dictionary, see `lada.jobs.JobDict`) is
    beyond the scope of this tutorial. 

    Prep-ing up: Exploring and modifying a job-dictionary before launch
    -------------------------------------------------------------------

    The following expects that you have a job-dictionary *saved to file* in
    some directory. How one gets a job-dictionary is somewhat dependent on how
    it is built, e.g. from a script, or from scratch using the capabilities
    described in `lada.jobs`. 

    >>> %explore path/to/jobdictionary # opens file *jobdictionary* in directory path/to/

    The command above opens a job-dictionary exactly as it exists on disk. It
    also keeps track of the path so that the job-dictionary can be save after
    being modified:

    >>> %savejobs                      # saves job-dictionary to the file it was read from.
    >>> %savejobs new/path/to/newfile  # saves to new file

    Note that in the first case above, the dictionary on disk will be
    overwritten. This may be good or bad... 
    
    Lets take up the example calculations given above. The job-dictionary would
    have the following tree-like structure:

    ::

       /           <-- root of the job dictionary
        ZnMgO/
          Spinel/  <-- parameters for calculation at this level
          Olivine/ <-- and here as well
        ZnRhO/
          Spinel/  <-- and here
          Olivine/ <-- and here

    When opening the job-dictionary, you are at the root level '/'. The list of
    all sub-jobs (eg subdirectories) can be printed using *%listjob*, the
    equivalent of *ls* for job-dictionaries:

    >>> %listjob
    ZnMgO ZnRhO

    To navigate to, say, ZnMgO, we can use *goto*, which is the equivalent of *cd*:

    >>> %goto ZnMgO
    >>> %listjob
    Spinel Olivine

    We can also navigate back, or forward two, and back to the root:
    
    >>> %goto ..            # goes one back, the same way "cd .." does.
    >>> %goto ZnMgO/Olivine # goes down two directories
    >>> %goto /             # goes back to root.
    >>> %listjob 
    ZnMgO ZnRhO

    Ok, so at this point you are not sure where the hell you are in the job
    dictionary. The equivalent of ``pwd`` is again ``goto`` (or ``explore``)
    but without any arguments:

    >>> %goto
    Current position in job dictionary: / 
    Filename of jobdictionary: path/to/jobdict

    Note that tab-completion works for all these commands. Try it!  Once we
    have navigated to a subjob where an actual calculation exists -- eg
    /ZnMgO/Olivine, as opposed to /ZnMgO -- we can edit both the crystal
    structure and the functional. 

    >>> %explore structure 
    >>> %explore functional 

    The first command opens an editor(vim? emacs? others?  see `Installing the
    lada interface to ipython.`_) with the POSCAR file corresponding to the
    structure. Once you save the file and exit the editor, the POSCAR is
    interpreted and replaces the old structure. Similarly, the second command
    will start an editor with python code:

    >>> from lada.vasp import Vasp
    >>> from lada.vasp.incar._params import NElect, UParams, Precision, Encut, Algo, Ediff
    >>> functional = Vasp()
    >>> functional.nelect = None
    >>> functional.encut = Encut(1)
    >>> ...

    >>> ...
    >>> # Parameters in the functional above **will** be
    >>> # overwritten by the following corresponding parameters.
    >>> jobparams = {}
    >>> jobparams["nelect"] = 1

    As shown above, there are two parts to this file. In the first part, the
    functional is defined. The name ``functional`` is crucial. It is simply the
    `vasp.Vasp` wrapper in this case, but could somewhat more complex, say a
    *meta*-functional over Vasp, such as `vasp.method.RelaxCellShape`, which
    would perform more than one calculation (eg relaxation steps + final static
    calculation). Whatever is called ``functional`` once the file interpreted
    by python will end up being the functional called by the job.  The second
    part are parameter specific to this particular calculation: parameters
    defined in jobparams will override parameters with the same name in the
    functional. In the example above, modifying ``nelect`` in the top part will
    have no effect (unless you remove it from the bottom). Modifying ``encut``
    in the top half will work however, since it is not present in the bottom
    half. This setup happens because of the way most job-dictionaries are
    constructed, using a template functional which is then modified to suit the
    particular purpose of a particular calculation. Note that this true only at
    construction time: modifying the functional when calling ``showme
    functional`` will **not** affect any other jobs.
    
    Launching highthoughput calculations
    ------------------------------------

    Once the job-dictionary has been modified, it can be launched with a single
    command:

    >>> # first load the dictionary
    >>> %explore path/to/dict
    >>> # then have fun with it
    >>> ...
    >>> # finally launch it.
    >>> %launch scattered 

    The last line will launch one `PBS
    <http://www.mcs.anl.gov/research/projects/openpbs/>`__ (or `SLURM
    <https://computing.llnl.gov/linux/slurm/>`__) job per actual calculation.
    That is the meaning of ``scattered``. Other ways of launching calculations
    will come as needed. The line above will allocate as many processors as
    there are atoms in the structure. This can be changed to suit you taste.
    Details on how to proceed are documented by ``launch`` itself.

    >>> # description of possible launch scenarios.
    >>> %launch --help

    >>> # description of possible *scattered* scenarios.
    >>> %launch scattered --help


    Inspecting and modifying jobs which did not complete:
    -----------------------------------------------------

    If all goes well, you won't need this. But life's a beautiful b*ch. At any
    time after launching a job-dictionary, you can inspect jobs that did not
    complete using the command:

    >>> %explore errors path/to/jobdictionary

    This will load the requested job-dictionary and mark those jobs wich are
    either running (on redmesa only) and jobs that completed successfully as
    uninteresting. You are left with errors only (note that the description of
    the other jobs is not lost and can be easily recovered if necessary). You
    can then rapidly navigate through the errors using the ``goto`` command we
    examined previously.

    >>> %goto next        # first failed job
    >>> vi OUTCAR         # diagnose
    >>> vi stderr         # diagnose
    >>> %showme functional # edit and correct.

    This will got the first error and *change the working directory* so you can
    inspect the output of the failed job. At this time, since the jobs has been
    launched, a directory on disk exists which reflects the "job"-directory of
    the job-dictionary.

    You can then go to the next failed job.
    
    >>> %goto next        # second failed job
    >>> vi OUTCAR         # diagnose
    >>> vi stderr         # diagnose
    >>> %showme functional # edit and correct.

    When no failed jobs are left, you will get the following message:

    >>> %goto next        
    Reached end of job list.

    If needed, you go through the failed jobs again:

    >>> %goto reset # restart the merry go round.

    In some cases, errors may happen sooner than we would like, and a printout
    only available in the standard error/ouput of the PBS job itself.

    >>> %goto pbs

    Will relocate you to the directory where all the pbs scripts and output are
    located.

    At this point, you can relaunch the failed jobs the same way you did
    previously.

    >>> %launch scattered

    Note that this will overwrite (don't worry, it will prompt you first) the
    job-dictionary on disk. Whether this is a problem is up to you. *None* of
    the information about any calculation has been lost, only information about
    the failed calculations may have been modified (by you). In any case, the
    modified job-dictionary can always be save to a different file with
    ``savejobs path/to/newdict``.


    Inspecting successful jobs and gathering output:
    ------------------------------------------------

    Successfull jobs are inspected the same way errors were above:

    >>> %explore results path/to/job

    Don't worry if you did overwrite the job-dictionary in the previous
    section. Didn't I say no job information was lost?  Jobs can be cycled
    through on a case by case basis using ``goto next``.  It is also possible
    to gather all data from a job, say eigenvalues

    >>> collect.eigenvalues
    { "ZnMgO/Spinel": array([[6.5,...,18.546513]]) * eV,
      ...
      "ZnRhO/Olivine": array([[0.1541631, ..., 0.1315]]) * eV }

    This results are arranged in  python `dict`_ where each key points to a
    particular job. Most quantities are numpy_ arrays_. `numpy`_ is the
    standard python package for scientific calculations. Furthermore, many
    quantities are given physical units using the quantities_ package.
    Gathering data may take a while *the first time around*, depending on the
    size of the job-dictionary. However, results are cached and will show up in
    a jiffy the second time around.

    .. _dict : http://docs.python.org/tutorial/datastructures.html#dictionaries
    .. _numpy : http://docs.scipy.org/doc/numpy/reference/
    .. _arrays : http://docs.scipy.org/doc/numpy/reference/arrays.html
    .. _quantities : http://packages.python.org/quantities/user/tutorial.html

    Which quantities are available, you ask? just hit tab at the right spot and
    look for available completions.

    >>> collect. # hit tab now. yes there is a dot at the end.

    Lots of printout, but if it doesn't sound like a physical property, use at
    your own risk.  Note that if data gathering encounters an error in any
    given job, it will fail silently for that job. If a particular calculation
    does not show up in the dictionary, either that quantity is not available
    or something went wrong in that job.

    Finally, when changing to particular leaf of the jod-dictionary using
    ``goto``, only subjobs will appear. When there is only one job, the
    quantities are returned explicitely, without the help of a dictionary.

    >>> %goto /ZnMgO
    >>> collect.eigenvalues
    { "ZnMgO/Spinel": array([[6.5,...,18.546513]]) * eV,
      "ZnMgO/Olivine": array([[-12.25463, ..., 21.216515]]) * eV }
    >>> %goto Spinel
    >>> collect.eigenvalues
    array([[6.5,...,18.546513]]) * eV


    Inspecting running job.
    -----------------------

    >>> %explore running path/to/job
    
    Well, that is there to. At least on redmesa.


    Example High-throughput script: point-defects.
    ----------------------------------------------

    The following describes, succintly, how to get high-throuput point-defects
    up and running.

    The typical hight-thoughput script is organized around to files: 
    
    - input.py : an input file where all relevant parameters are gathered.
    - test.py  : the script itself. 

    Here is how you can retrieve and set up these scripts on **redmesa**:

    >>> mkdir ladatest
    >>> cd ladatest
    >>> cp ~mmdavezac/usr/src/LaDa/master/test/point_defects/*.py .

    Since the high-thoughput for point-defects uses Vasp, we will also need
    pseudo-potential files as input. Those for this specific test can be found
    here: 

    >>> cp -r ~mmdavezac/usr/src/LaDa/master/test/point_defects/pseudos .

    Now that we all the files, we can start the IPython interface:

    >>> ~mmdavezac/bin/pyladas.sh

    This should launch IPython with all paths correctly initialized.
    At this point, you could edit the ``input.py`` file to suit your needs. If
    endowed with bravery, you are highly encourage to edit the ``test.py``
    script itself and play with it. 

    The ``input.py`` file is all python code. It will be interpretated by
    python when running ``test.py``. Anything python can be done in this file.
    All that matters is that at the end of the day a few variables have been
    defined. Lets look at a few, in order of appearance.

    >>>  from lada.crystal import A2BX4
    >>>  lattice = A2BX4.b5()

    >>>  for site in lattice.sites:
    >>>    site.type = {"A": "Rh", "B": "Zn", "X": "O"}[site.type[0]]
    >>>  lattice.scale = 8.506

    These two lines create a spinel lattice (b5) using the stock A2BX4 lattices
    from lada. The second set of lines replaces the stock names ``A``, ``B``,
    ``X`` with appropriate elements and sets the lattice parameter.
    
    
    You could replace this with any other stock lattice in `lada.crystal.A2BX4`
    or `lada.crystal.bravais`, or one of your own making (see
    `lada.crystal.A2BX4` source code for templates on how to do this). You
    could also read a poscar using `lada.crystal.read_poscar`, and then
    transform the returned variable to a lattice using: 

    >>> from lada.crystal import read_poscar, structure_to_lattice
    >>> structure = read_poscar(["Rh", "Zn", "O"], "path/to/poscar")
    >>> lattice = structure_to_lattice(structure)
    

    The second variable to be defined is the supercell to use for
    point-defects.  This is a simple `numpy`_ array. 

    >>>  supercell = array([[1, 0, 0],
    >>>                     [0, 1, 0],
    >>>                     [0, 0, 1]], dtype="float64" )
       

    The third type of variable is the vasp functional

    >>>  vasp = Vasp()
    >>>  vasp.kpoints    = "Automatic generation\\n0\\nMonkhorst\\n2 2 2\\n0 0 0"
    >>>  vasp.precision  = "accurate"
    >>>  vasp.ediff      = 1e-5
    >>>  vasp.encut      = 1
    >>>  ...

    Don't forget to include all relevant pseudo-potentials!
      
    The fourth variable is somewhat more bizare. It is python `dict`_ where
    each key/value pair will override a vasp parameter during the
    ionic-relaxation step. Indeed, each point-defect structure if fist relaxed,
    and then recomputed statically for maximum accuracy. Be default, the
    dictionary is empty: nothing is overriden.

    >>>  relaxation_parameters = {}
      

    Finally getting to the meat of it. Substitutional and vacancy defects are
    defined in the following line:

    >>>  substitutions = {"Rh": ["Zn", None], "Zn": ["Rh", None], "O": [None]}

    This reads as in: compute Zn on Rh, vacancy on Rh, Rh on Zn, vacancy on Zn,
    vacancy on O.
      
    The interstitials are a bit more complicated to define, since an actual
    position in cartesian coordinates must be given, as well as a name.

    >>>  interstitials = { "Li": [(0,0,0, "16c"), (0.75,0.75,0.75, "32e_0.75")],
    >>>                    "Na": [(0,0,0, "16c"), (0.75,0.75,0.75, "32e_0.75")] }

    It is important that the name of the interstial defect differ for any given elements. 



    Once the input has been modified to suit your needs, you can create the
    jobdictionary (within pylada) with 

    >>> import test
    >>> jobdict = test.create_jobdict("input.py")
    >>> %explore jobdict
    >>> %savejobs path/to/file/on/disk

    - The first line imports the test script in memory.
    - The second line creates the job-dictionary in memory.
    - The third line makes ``jobdict`` the current job-dictionary, still only in
      memory (The current job-dictionary is ``current_jobdict`` in your IPython
      user namespace).
    - The fourth line save the current job-dictionary to file. 

    At this point, you can examine and launch the job-dictionary, as explainged
    in sections `Prep-ing up: Exploring and modifying a job-dictionary before
    launch`_ and `Launching highthoughput calculations`_.



"""
__docformat__ = "restructuredtext en"
from contextlib  import contextmanager
from _collect import Collect

def _get_current_job_params(self, verbose=0):
  """ Returns a tuple with current job, filename, directory. """

  ip = self.api
  if "current_jobdict" not in ip.user_ns: 
    if verbose > 0:
      print "No job dictionary. Please use \"explore\" magic function."
    return None, None
  current = ip.user_ns["current_jobdict"]
  if "current_jobdict_path" not in ip.user_ns: 
    if verbose > 1:
      print "No filepath for current job dictionary.\n"\
            "Please set current_jobdict_path."
    return current, None
  path = ip.user_ns["current_jobdict_path"]
  return current, path

def listjobs(self, arg):
  """ Lists subjobs. """
  ip = self.api
  current, path = _get_current_job_params(self, 1)
  ip.user_ns.pop("_lada_error", None)
  if current == None: return
  if len(arg) != 0:
    if arg == "all": 
      for job, d in current.root.walk_through():
        if job.is_tagged: continue
        print job.name
      return
    try: subdict = current[arg] 
    except KeyError:
      print "%s is not a valid jobname of current job dictionary." % (arg)
      return
    current = current[arg]
    if not hasattr(current, "children"):  
      print "%s is not a valid jobname of current job dictionary." % (arg)
      return
  if len(current.children) == 0: return
  string = ""
  lines = ""
  for j in current.children.keys():
    if current.children[j].is_tagged: continue
    if len(string+j) > 60:
      if len(lines) != 0: lines += "\n" + string
      else: lines = string
      string = ""
    string += j + " "

  if len(lines) == 0: print string
  else: print lines + "\n" + string


def saveto(self, event):
  """ Saves current job to current filename and directory. """
  from os.path import exists, abspath, isfile
  from lada import jobs
  ip = self.api
  # gets dictionary, path.
  current, path = _get_current_job_params(self, 1)
  ip.user_ns.pop("_lada_error", None)
  if current == None:
    ip.user_ns["_lada_error"] = "No job-dictionary to save."
    print ip.user_ns["_lada_error"] 
    return
  args = [u for u in event.split() ]
  if len(args) == 0: 
    if path == None: 
      ip.user_ns["_lada_error"] = "No current job-dictionary path.\n"\
                                  "Please specify on input, eg"\
                                  ">saveto this/path/filename"
      print ip.user_ns["_lada_error"] 
      return
    if exists(path): 
      if not isfile(path): 
        ip.user_ns["_lada_error"] = "%s is not a file." % (path)
        print ip.user_ns["_lada_error"] 
        return
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("File %s already exists.\nOverwrite? [y/n] " % (path))
      if a == 'n':
       ip.user_ns["_lada_error"] = "User said no save."
       return
    jobs.save(current.root, path, overwrite=True) 
  elif len(args) == 1:
    if exists(args[0]): 
      if not isfile(args[0]): 
        ip.user_ns["_lada_error"] = "%s is not a file." % (path)
        print ip.user_ns["_lada_error"] 
        return
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("File %s already exists.\nOverwrite? [y/n] " % (args[0]))
      if a == 'n':
       ip.user_ns["_lada_error"] = "User said no save."
       return
    jobs.save(current.root, args[0], overwrite=True) 
    ip.user_ns["current_jobdict_path"] = abspath(args[0])
    if "collect" not in ip.user_ns: ip.user_ns["collect"] = Collect()
  else:
    ip.user_ns["_lada_error"] = "Invalid call to saveto."
    print ip.user_ns["_lada_error"] 


def current_jobname(self, arg):
  """ Returns current jobname. """
  ip = self.api
  if "current_jobdict" not in ip.user_ns: return
  print ip.user_ns["current_jobdict"].name
  return

def fakerun(self, event):
  """ Creates job directory tree and input files without computing. """
  from os.path import split as splitpath, exists, isdir
  ip = self.api

  current, path = _get_current_job_params(self, 2)
  ip.user_ns.pop("_lada_error", None)
  if current == None or path == None: return
  if len(event.split()) > 1: 
    print "fakerun does not take an argument."
    return
  elif len(event.split()) == 1: directory = event.split()[0]
  else: directory = splitpath(path)[0]

  if exists(directory) and not isdir(directory):
    print "%s exists and is not a directory." % (directory)
    return 
  elif exists(directory):
    a = ''
    while a not in ['n', 'y']:
      a = raw_input("%s exists. \n"\
                    "Some input files could be overwritten.\n"\
                    "Continue? [y/n]" % (directory))
    if a == 'n': return
  for job, dirname in current.walk_through(directory):
    if not job.is_tagged: job.compute(outdir=dirname, norun=True)

def run_current_jobdict(self, event):
  """ Runs job dictionary interactively. """
  from os.path import split as splitpath, exists, isdir
  ip = self.api

  current, path = _get_current_job_params(self, 2)
  ip.user_ns.pop("_lada_error", None)
  if current == None or path == None: return
  if len(event.split()) > 1: 
    print "fakerun does not take an argument."
    return
  elif len(event.split()) == 1: directory = event.split()[0]
  else: directory = splitpath(path)[0]

  if exists(directory) and not isdir(directory):
    print "%s exists and is not a directory." % (directory)
    return 
  elif exists(directory):
    a = ''
    while a not in ['n', 'y']:
      a = raw_input("%s exists. \n"\
                    "Some input files could be overwritten.\n"\
                    "Continue? [y/n]" % (directory))
    if a == 'n': return
  for job, dirname in current.walk_through(directory):
    if not job.is_tagged: job.compute(outdir=dirname)

def qstat(self, arg):
  """ squeue --user=`whoami` -o "%7i %.3C %3t  --   %50j" """
  from subprocess import Popen, PIPE
  from IPython.genutils import SList

  ip = self.api
  # finds user name.
  whoami = Popen(["whoami"], stdout=PIPE).stdout.readline()[:-1]
  squeue = Popen(["squeue", "--user=" + whoami, "-o", "\"%7i %.3C %3t    %j\""],
                 stdout=PIPE)
  result = squeue.stdout.read().rstrip().split('\n')
  result = SList([u[1:-1] for u in result])
  return result.grep(str(arg[1:-1]))

def cancel_completer(self, info):
  return qstat(self, info.symbol).fields(-1)[1:]

def cancel_jobs(self, arg):
  """ Cancel jobs which grep for whatever is in arg.
  
      For instance, the following cancels all jobs with "anti-ferro" in their
      name.

      >>> %cancel_jobs "anti-ferro"
  """
  from subprocess import Popen, PIPE
  
  arg = str(arg[1:-1])
  if len(arg) == 0: 
    print "cancel_job Requires an argument."
    print "Please use please_cancel_all_jobs to cancel all jobs."
    return
  result = qstat(self, arg)
  for u, name in zip(result.fields(0), result.fields(-1)):
    print "cancelling %s." % (name)
  a = ''
  while a not in ['n', 'y']:
    a = raw_input("Are you sure you want to cancel the jobs listed above? [y/n] ")
  if a == 'n': return
  for u, name in zip(result.fields(0), result.fields(-1)):
    self.api.system("scancel %i" % (int(u)))

def please_cancel_all_jobs(self, arg):
  """ Cancel all jobs. """
  from subprocess import Popen, PIPE
  
  a = ''
  while a not in ['n', 'y']: a = raw_input("Are you sure you want to cancel all jobs? [y/n] ")
  if a == 'n': return
  result = qstat(self, None)
  for u in result.field(0):
    self.api.system("scancel %i" % (int(u)))

def ipy_init():
  """ Initialises ipython session. 

      In order to initialize an ipython session with the lada interface, the
      following lines can be introduces in the main function of your
      ipy_user_conf.py

      >>> try: import lada.ipython 
      >>> except ImportError:
      >>>   print "Could not find module lada.ipython."
      >>>   pass
      >>> else: lada.ipython.ipy_init()
  """ 
  try: import IPython.ipapi
  except: pass
  else:
    from os import environ
    import lada
    from ._goto import goto, iterate, goto_completer
    from ._explore import explore, explore_completer
    from ._showme import showme, showme_completer
    from ._launch import launch, launch_completer
    
    ip = IPython.ipapi.get()
    ip.expose_magic("explore", explore)
    ip.expose_magic("goto", goto)
    ip.expose_magic("listjobs", listjobs)
    ip.expose_magic("jobname", current_jobname)
    ip.expose_magic("iterate", iterate)
    ip.expose_magic("showme", showme)
    ip.expose_magic("savejobs", saveto)
    ip.expose_magic("fakerun", fakerun)
    ip.expose_magic("launch", launch)
    ip.expose_magic("run_current_jobdict", run_current_jobdict)
    ip.set_hook('complete_command', goto_completer, re_key = '\s*%?goto')
    ip.set_hook('complete_command', showme_completer, re_key = '\s*%?showme')
    ip.set_hook('complete_command', explore_completer, re_key = '\s*%?explore')
    ip.set_hook('complete_command', launch_completer, re_key = '\s*%?launch')
    if "SNLCLUSTER" in environ:
      if environ["SNLCLUSTER"] in ["redrock", "redmesa"]:
        ip.expose_magic("qstat", qstat)
        ip.expose_magic("cancel_jobs", cancel_jobs)
        ip.set_hook('complete_command', cancel_completer, re_key = '\s*%?cancel_jobs')
        ip.expose_magic("please_cancel_all_jobs", please_cancel_all_jobs)
    
    for key in lada.__all__:
      if key[0] == '_': continue
      if key == "ipython": continue
      if key == "jobs": ip.ex("from lada import jobs as ladajobs")
      else: ip.ex("from lada import " + key)

