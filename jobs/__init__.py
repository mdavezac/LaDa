""" Allows the creation of jobs. 

    Contains the following methods and classes.
       - JobDict: tree of calculations and subcalculations.
       - walk_through: iterates over calculation in a tree.
       - save: pickles a tree to file. Acquires a file lock to do this! If you
           always use a file locking mechanism, then this should work out of the
           box for you.
       - load: loads a pickled tree from file. Acquires a file lock to do this! If you
           always use a file locking mechanism, then this should work out of the
           box for you.
       - pbs_script: Creates pbs-script to perform calculations on a tree.
"""
from ..opt.decorators import add_setter, broadcast_result

class JobDict(object):
  """ Tree of jobs. 

      A tree of jobs can be created as:

        >>> jobtree = JobDict() # root tree
        >>> subjob1 = jobtree / "subjobname1"
        >>> subjob2 = jobtree / "subjobname2"

      This creates a two sub-jobs called subjobname1 and subjobname2 (bothe are strings.).
      Each subjob can be further expanded as:

        >>> subsubjob = subjob1 / "subsubjobname"
        >>> subsubjob = jobtree / "subjobname1" / "subsubjobname"

      At this point however, these jobs do not do anything. They haven't been
      given anything to do:

        >>> subjob1.functional = some callable object.

      Where a callable object is anything which can be called (eg, as in
      C{method(parameters)}), C{method} is the callable object:

        >>> subjob1.functional = method

      Note that if you want to pickle the jobtree (eg save to file, google
      "python pickle"), then method should be pickleable. 
      The argument to the method (if any) are specified as:

        >>> subjob1.args = (arg1) # one argument. Not the parenthesis. This is a tuple.
        >>> subjob1.args = (arg1, arg2) # two arguments.

      Keyword arguments can be specified as:

        >>> subjob1.jobparams["keywordname"] = value

      Once they have been specified, keyword arguments can be accessed directly:

        >>> # first specify
        >>> subjob1.jobparams["keywordname"] = value1
        >>> # then change
        >>> subjob1.keywordname = value2

      Where keywordname in the last line should be substituted for its actual
      value (eg structure, nelect...)

      Having a jobtree is great, looping through it is better:

        >>> for job, name in jobtree.walkthrough("root_result_directory"):
        >>>    result = job.compute(outdir=name)


      The L{walkthrough} method does just that: it goes through the whole tree
      returning those branches where there is something to execute (eg
      C{job.functional != None}). It also returns outdir which is the subdirectory
      where to execute that particular job. Using the jobtree defined above,
      for which only subjob1 was actually given something to do, 

        >>> for job, name in walkthrough(jobdict=jobtree, outdir="root_result_directory"):
        >>>    print name

      would print out: "root_result_directory/subjoname1".
      The actual directory is not created. The name is only a suggestion to use or not.
        
      Putting it all together:

        >>> # create a job tree
        >>> jobtree = JobDict() # root tree
        >>> subjob1 = jobtree / "subjobname1"
        >>> subjob2 = jobtree / "subjobname2"
        >>> subsubjob = jobtree / "subjobname1" / "subsubjobname"
        >>> # specify an actual job in one of the subjobs.
        >>> subjob1.functional = method
        >>> subjob1.args = (arg1, arg2) # two arguments.
        >>> subjob1.add_param = "whatever", value1
        >>> subjob1.whatever = value2
        >>> # actually performs calculations.
        >>> for job, name in walkthrough(jobdict=jobtree, outdir="root_result_directory"):
        >>>    job.compute(outdir=name)

      Using the jobtree defined above, this would end-up calling:

        >>> metod(arg1, arg2, whatever=value1, outdir="root_result_directory/subjoname1") 


      The code above does not work if more than one pbs script is launched. In that case, we want:

        >>> # somewhere the jobdict was created.
        >>> jobdict.save("pickle")
        >>> # somewhere else, most likely in another script, the jobs are executed.
        >>> for job, name in jobs.bleed(comm=local_comm):
        >>>    result = job.compute(outdir=name)

      L{bleed} makes sure that once a job is accepted, it will not be accessed
      by any other process or pools of processes. L{bleed} allows parallelization
      over pbs scripts and pools of MPI processes. However, it does change the
      C{jobdict} saved to "pickle". This means that once all jobs are executed,
      L{bleed} (but not L{walk_through}) will find that C{jobdict} is empty. To
      undo these changes and use the jobdict for, say, ouptut analysis, simply
      use L{walk_through} where possible (single pbs script or interactive
      sessions), or do:

      >>> jobs.unbleed("pickle")

      Beware, at this point, relaunching the job could overwrite, if the
      functional allows it... In practice, bled and unbled jobs can also be
      selected with:

      >>> for job, outdir in job.walk_trough():
      >>>   if job.is_tagged: # this job has been bled/tagged.
      >>>        # do something.
      >>>        job.untag() # this job will no longuer be considered bled/tagged.
      >>>   else:  # this job was not previously bled/tagged.
      >>>        # do something.
      >>>        job.tag() # this job will now be considered bled/tagged.
     
      Coding: JobDict has the following attributes:
        - children: A dict object holding instances of JobDict. These are the sub-jobs.
        - jobparams: All parameters regarding actual calculations. It contains,
          at start, only two predefined parameters.
           - functional: is the callable (preferably pickleable) to execute.
           - args: is a tuple of arguments (eg. C{functional(*args)}
           - all others are keyword arguments.
        - parent is an instance to the parent job (eg the instance which holds
          self in children) or None.
        - It may also have a _tagged attribute to check for bled/unbled jobs.
      The __getattr__, __setattr__, and __delattr__ have been rewired to
      perform on objects in jobparams. Note however that __setattr__ will not
      set new object in jobparams, but rather pass on the call to the parent
      class' __setattr__. To set new job parameters, one should use
      L{add_param} or jobparams directly.
  """

  def __init__(self):
    super(JobDict, self).__init__()
    # List of subjobs (as in subdirectories). 
    super(JobDict, self).__setattr__("children", {})
    # This particular job. 
    super(JobDict, self).__setattr__("jobparams", {})
    # Parent job. 
    super(JobDict, self).__setattr__("parent", None)

    # no jobs yet.
    self.jobparams["functional"] = None
    # no arguments yet.
    self.jobparams["args"] = (None)
    # no restart parameters yet.
    self.jobparams["restart"] = None
    
  @property
  def name(self):
     """ Returns the name of this dictionary, relative to root. """
     if self.parent == None: return "/"
     string = None
     for key, item in self.parent.children.items():
       if id(item) == id(self):
         string = self.parent.name + key
         break
     assert string != None, RuntimeError("Could not determine the name of the dictionary.")
     if not self.is_job: string += "/"
     return string

  def __getitem__(self, index): 
    """ Returns job description from the dictionary.

        If the job does not exist, will create it.
    """
    from re import split
    from os.path import normpath

    index = normpath(index)
    if index == "" or index == None or index == ".": return self
    if index[0] == "/": return self.root[index[1:]]

    result = self
    names = split(r"(?<!\\)/", index)
    for i, name in enumerate(names):
      if name == "..":
        if result.parent == None: raise KeyError("Cannot go below root level.")
        result = result.parent
      elif name in result.children: result = result.children[name]
      elif name in result.jobparams:
        if i+1 != len(names): raise KeyError("job or job parameter " + index + " does not exist.") 
        return result.jobparams[name]
      else: raise KeyError("job or job parameter " + index + " does not exist.")
    return result
 
  def __delitem__(self, index): 
    """ Returns job description from the dictionary.

        If the job does not exist, will create it.
    """
    from os.path import normpath, relpath

    index = normpath(index)

    try: deletee = self.__getitem__(index) # checks if exists.
    except KeyError: raise

    if isinstance(deletee, JobDict): 
      assert id(self) != id(deletee), KeyError("Will not commit suicide.")
      parent = self.parent
      while parent != None: 
        assert id(parent) != id(deletee), KeyError("Will not go Oedipus on you.")
        parent = parent.parent

    parent = self[index+"/.."]
    name = relpath(index, index+"/..")
    if name in parent.children:
      assert id(self) != id(parent.children[name]),\
             KeyError("Will not delete self.")
      return parent.children.pop(name)
    if name in parent.jobparams: return parent.jobparams.pop(name)
    raise KeyError("job or job parameter " + index + " does not exist.")

  def __setitem__(self, index, value): 
    """ Sets job/subjob description in the dictionary.
    
        If the job does not exist, will create it.
    """
    from re import split
    from os.path import normpath, relpath

    index = normpath(index)
    assert index not in ["", ".", None], KeyError("Will not set self.")
    assert index[0], KeyError("Will not set root: " + index + ".")

    result = self.__div__(index+"/..")
    name = relpath(index, index+"/..")
    if name in result.children  and isinstance(value, JobDict): 
                                     result.children [name] = value
    elif name in result.jobparams:   result.jobparams[name] = value
    elif isinstance(value, JobDict): result.children [name] = value
    else:                            result.jobparams[name] = value

  def __div__(self, index): 
    """ Adds name as a subtree of self. """
    from re import split
    from os.path import normpath

    index = normpath(index)
    if index in ["", ".", None]: return self
    if index[0] == "/":  # could create infinit loop.
      result = self
      while result.parent != None: result = result.parent
      return result / index[1:]

    names = split(r"(?<!\\)/", index) 
    result = self
    for name in names:
      if name == "..":
        if result.parent != None: result = result.parent
        continue
      elif name not in result.children:
        result.children[name] = JobDict()
        result.children[name].parent = result
      result = result.children[name]
    return result

  @property
  def is_job(self):
    """ True if self.job has keyword \"functional\". """
    if "functional" not in self.jobparams: return False
    return self.jobparams["functional"] != None

  def items(self):
    """ Iterator over children jobs. """
    return self.children.items()
  def subjobs(self):
    """ Iterator over children jobs. """
    return sorted(self.children.keys())
    
  def compute(self, **kwargs):
    """ Performs calculations over job list. """  

    kwargs.update(self.jobparams)
    if "functional" not in kwargs: return None
    functional = kwargs.pop("functional")
    if functional == None: return
    args = kwargs.pop("args", ())
    # restart is bit painful to deal with, so farm it out to private function for now.
    if "restart" in kwargs:
      args, kwargs = self._restart(*args, **kwargs)

    if args == None: args = ()
    assert hasattr(args, "__iter__"),\
           RuntimeError("Functional argument \"args\" is not a sequence.")
    return functional(*args, **kwargs)

  def _restart(self, *args, **kwargs):
    """ Restart for functionals.
    
        Transforms arguments and dictionary to deal with restart request.
        In this case, modifies the first argument of the functional. Then  sets
        keyword argument restart to None (so nothing else is restarted).
    """
    from lada.opt.changedir import Changedir
    from sys import stderr
    if kwargs["restart"] == None: return args, kwargs

    restart = kwargs.pop("restart")
    kwargs["restart"] = None

    outdir = kwargs["outdir"] if "outdir" in kwargs else "."
    comm = kwargs["comm"] if "comm" in kwargs else None
    with Changedir(outdir) as pwd:
      print "restart: ", outdir, restart[1]
      restart = restart[0](restart[1], comm=comm)
    if not restart.success:
      # cannot perform this job since dependency is not successfull
      comm = kwargs["comm"] if "comm" in kwargs else None
      if (comm.rank == 0 if comm != None else True):
        print >> stderr, "Could not perform job %s since it depends upon %s completing first."\
              % (outdir, join(outdir, restart[1]))
      class NoSuccess:
        def __init__(self): self.success, self.directory = False, outdir
      return NoSuccess()
    args[0] = restart.structure
    return args, kwargs

  def update(self, other):
    """ Updates job and tree with other.
    
        Updates the dictionaries of job parameters and sub-jobs. Will overwrite
        if items in C{other} are found in self.
        @param other: An other job tree from which to update.
        @type other: L{JobDict}
    """
    self.children.update(other.children)
    self.jobparams.update(other.jobparams)

  def __str__(self):
    result = "Jobs: \n"
    for dummy, name in self.walk_through():
      result += "  " + name + "\n"
    return result

  @property
  def untagged_jobs(self):
    """ Returns a string with only untagged jobs. """
    result = "Jobs: \n"
    for dummy, name in self.walk_through():
      if not dummy.is_tagged: 
        result += "  " + name + "\n"
    return result

  @property
  def is_tagged(self): return hasattr(self, "_tagged")

  def tag(self):
    """ Adds _tagged attribute. """
    if self.is_job: super(JobDict, self).__setattr__("_tagged", True)
    
  def untag(self):
    """ Adds _tagged attribute. """
    if hasattr(self, "_tagged"): self.__delattr__("_tagged")

  def __delattr__(self, name):
    """ Deletes job attribute. """
    if name in self.jobparams: return self.jobparams.pop(name)
    elif name in self.__dict__: return self.__dict__.pop(name)
    raise AttributeError("Unknown job attribute " + name + ".")

  def __getattr__(self, name):
    """ Returns job attribute. """
    if name in self.jobparams: return self.jobparams[name]
    raise AttributeError("Unknown job attribute " + name + ".")

  def __setattr__(self, name, value):
    """ Sets job attribute. """
    if name in self.jobparams: self.jobparams[name] = value
    else: super(JobDict, self).__setattr__(name, value)

  @add_setter
  def add_param(self, args):
    """ Adds a job parameter.
    
        Actual job parameter must be first set as:

         >>> jobs.add_param = "name", value

        Then they can be accessed directly as:

         >>> jobs.name = othervalue

        Where name should be the value of the string given in the first line.
        The only preset job parameter are "functional" and "args", as stand in for
        the functional to use and the tuple to use as its arguments.
    """
    assert len(args) < 3, ValueError("Too many parameters given to add_param: %s." % (args))
    name = args[0]
    value = None if len(args) == 1 else args[1]
    self.jobparams[name] = value

  def __getstate__(self):
    d = self.__dict__.copy()
    params = d.pop("jobparams")
    return d, params, "This is a JobDict pickle. Grep me!"
  def __setstate__(self, args):
    super(JobDict, self).__setattr__("jobparams", args[1])
    d = self.__dict__.update(args[0])

  def walk_through(self, outdir=None):
    """ Iterates over jobs. 

        @return: yields (job, directory):
          - job is a Jobdict which contains something to execute.
          - directory is a suggested directory name with C{outdir} as its root.
    """
    from os.path import join
    if outdir == None: outdir = ""
    # Yield this job if it exists.
    if self.is_job: yield self, outdir
    # Walk throught children jobdict.
    for name in self.subjobs():
      for u in self[name].walk_through(join(outdir, name)): 
        yield u

  @property
  def nbjobs(self):
    """ Returns the number of jobs in tree. """
    return len([0 for j, o in self.walk_through() if not j.is_tagged])

  @property 
  def root(self): 
    """ Returns root dictionary. """
    result = self
    while result.parent != None: result = result.parent
    return result

def walk_through(jobdict, outdir = None, comm = None):
  """ Generator to iterate over actual calculations. 
      
      see L{JobDict} description.
  """
  if isinstance(jobdict, str): jobdict = load(jobdict, comm=comm)
  for u in jobdict.walk_through(outdir): yield u

@broadcast_result(key=True)
def save(jobdict, path = None, overwrite=False, comm=None): 
  """ Pickles a job to file.

      This method first acquire an exclusive lock (using os dependent lockf) on
      the file before writing. This way not two processes can read/write to
      this file while using this function.
      @param jobdict: A jobtree to pickle. 
      @type jobdict: JobDict
      @param path: filename of file to which to save pickle. overwritten. If
        None then saves to "pickled_jobdict"
      @param comm: Only root process gets to do anything.
      @type comm: boost.mpi.communicator
      @param overwrite: if True, then overwrites file.
  """ 
  from os.path import exists
  from pickle import dump
  from ..opt import open_exclusive
  if path == None: path = "pickled_jobdict"
  if exists(path) and not overwrite: 
    print path, "exists. Please delete first if you want to save the job dictionary."
    return
  with open_exclusive(path, "wb") as file: dump(jobdict, file)
  print "Saved job dictionary to %s." % (path)

@broadcast_result(key=True)
def load(path = None, comm = None): 
  """ Unpickles a job from file.

      This method first acquire an exclusive lock (using os dependent lockf) on
      the file before reading. This way not two processes can read/write to
      this file while using this function.
      @param path: filename from which to load pickle. 
        If None then saves to "pickled_jobdict"
      @param comm: Broadcasts from root process. 
      @type comm: boost.mpi.communicator
      @return: Returns a JobDict object.
  """ 
  from os.path import exists
  from pickle import load as load_pickle
  from ..opt import open_exclusive
  if path == None: path = "pickled_jobdict"
  assert exists(path), IOError("File " + path + " does not exist.")
  with open_exclusive(path, "rb") as file: result = load_pickle(file)
  print "Loaded job list from", path, "."
  return result

def bleed(path=None, outdir=None, comm=None): 
  """ Generator which deepletes a job dictionary of its jobs. 

      This function alters the dictionary stored in C{path}. If C{path} is
      empty, then returns None, None. An exclusive lock is acquired before
      reading/writing to C{path}.  This way, if using L{bleed}, L{save},
      L{load}, two processes will not step on each others jobs.

      This function is different from JobDict.bleed in that it alters a
      dictionary stored in a file.
      @param: Filename of a pickled jobdictionary.
      @outdir: Root result directory. 
      @comm: Will broadcast yielded stuff from root. Because of file locking,
             this generator may freeze the system if not used correctly with mpi.
      @return: yields (job, directory), see L{walk_through}.
           - job: a job dictionary with the current job to execute.
           - directory: a suggested directory name with L{outdir} as its root.
  """
  from os.path import join, exists
  from pickle import load as load_pickle, dump
  from boost.mpi import broadcast
  from ..opt import acquire_lock
  if path == None: path = "pickled_jobdict"
  if outdir == None: outdir = ""

  is_root = True if comm == None else comm.rank == 0
  if is_root:
    while True:
      if not exists(path): 
        print "Job dictionary", path, "does not exist."
        return
      # acquires a lock file. 
      with acquire_lock(path) as lock_file:
        # tries to load pickle.
        with open(path, "rb") as file:
          try: jobdict = load_pickle(file)
          except EOFError: break
        # Checks if there are any jobs.
        if jobdict.nbjobs == 0: break
        # Pops first job.
        for job, directory in jobdict.walk_through():
          if not job.is_tagged: job.tag(); break
        # writes modified dictionary to path.
        with open(path, "wb") as file: dump(jobdict, file)
      broadcast(comm, (job, join(outdir, directory)), 0)
      yield job, join(outdir, directory)
    broadcast(comm, None, 0)
    return
  else: 
    while True:
      result =  broadcast(comm, root=0)
      if result == None: return
      else: yield result
  

def unbleed(path=None, comm=None): 
  """ Unbleeds a dictionary stored in a file. """
  from os.path import join, exists
  from pickle import load as load_pickle, dump
  from ..opt import acquire_lock

  # only one process needs do anything.
  if (True if comm==None else comm.rank == 0): 
    assert exists(path), IOError( "Job dictionary" + path + "does not exist.")
    with acquire_lock(path) as lock_file:
      with open(path, "rb") as file: jobdict = load_pickle(file)
      for job, outdir in jobdict.walk_through(): jobdict.untag()
      with open(path, "wb") as file: dump(jobdict, file)
  # synchronizes communicator.
  if comm != None: comm.barrier()


def unsucessfull(jobdict, extractor, outdir = None):
  """ Returns jobdictionary with unsucessfull/incomplete jobs.

      This function will not modify the input dictionary.
      @param jobdict: a job dictionary instance or the path to a pickled
        job-dictionary.
      @type jobdict: JobDict or str
      @param extractor: Some instance with a directory attribute and a success
        attribute, capable of judging the success of an operation.
      @return: A JobDict instance with incomplete and unsuccessful jobs.
  """
  from os.path import split as splitpath
  from copy import deepcopy

  extractor = deepcopy(extractor)

  # if path, get pickled dictionary.
  if not hasattr(jobdict, "unbleed"):
    if outdir == None and len(splitpath(jobdict)[0]) != 0:
      outdir = splitpath(jobdict)[0]
    jobdict = load(jobdict)
  # otherwise make deep copy.
  else: jobdict = deepcopy(jobdict)
  if outdir == None: outdir = ""

  # go through dictionary.
  for job, directory in jobdict.walk_through(outdir=outdir):
    extractor.directory = directory
    # and bleed successfull jobs.
    if extractor.success: job.tag()
    else: job.untag()

  return jobdict



def pbs_script( jobdict, outdir = None, template = None, \
                name = None, pickle = None, pyscript = None, **kwargs):
  """ Parallelizes jobs over different pbs scrits. 

      The root directory (outdir) is created if it does not exist, and the
      pbs-script placed there. The pbs-script is not launched. 

      @param jobdict: job-tree instance over which to parallelize.
      @param outdir: root directory of calculation. Current working directory if None.
      @type jobdict: L{JobDict}
      @param template: PBS-script template to use. See L{jobs.templates}.
        Default: L{jobs.templates.default_pbs}.
      @type template: callable.
      @param pickle: Filename to which L{JobDict} will be pickled.
      @param pyscript: Filename of python script to execute. Copied to outdir.
      @param kwargs: Passed on to template.
  """
  from shutil import copy
  from os import getcwd, makedirs, environ
  from os.path import abspath, join, exists, relpath, samefile, split as pathsplit
  from ..opt.changedir import Changedir
  from templates import default_pbs, default_slurm
 
  # sets up default input.
  if template == None:
    which = "SNLCLUSTER" in environ
    if which: which = environ["SNLCLUSTER"] in ["redrock", "redmesa"]
    template = default_slurm if which else default_pbs
  if outdir == None: outdir = getcwd() 
  if pickle == None: pickle = "job_pickle"
  if name == None: name = relpath(outdir, outdir+"/..")
  if pyscript == None: # just copy standard script.
    pyscript = __file__.replace(pathsplit(__file__)[1], "runme.py")
    if pyscript[-3:] == "pyc": pyscript = pyscript[:-1]
  assert exists(pyscript), RuntimeError("Execution script %s does not exist." % (pyscript))

  assert not exists(join(outdir, pickle)),\
         RuntimeError( "Job list %s already exist. Will not overwrite.\n"\
                       "Please stand back while I reboot the universe.\n"\
                       % (join(outdir, pickle)))


  # creates result dictionary.
  if not exists(outdir): makedirs(outdir)
  # pickle jobs.
  save(jobdict = jobdict, path = join(outdir, pickle))
  # copies script
  pyscript_filename = pathsplit(pyscript)[1]
  if not exists(join(outdir, pyscript_filename)): copy(pyscript, outdir)
  elif not samefile(pyscript, join(outdir, pyscript_filename)): copy(pyscript, outdir)
  # Writes out pbs scripts.
  result = []
  # writes pbs script.
  jobname = name
  with open(abspath(join(outdir, "launchme")), "w") as file:
    template( file, pickle=pickle, outdir=abspath(outdir),\
              name=jobname, pyscript=pyscript_filename, **kwargs)


def one_per_job(jobdict, outdir = None, mppalloc=None, ppath=None, **kwargs):
  """ Launches one pbs job per job. 
  
      @param jobdict: job dictionary.
      @param outdir: root output directory.
      @param mppalloc: an mpi allocation scheme. It takes the job as argument.
                       If a number, then flat allocation scheme across all jobs.
  """
  from os import environ, getcwd
  from os.path import abspath, join, split as pathsplit, expanduser
  from lada.opt.changedir import Changedir
  from templates import default_pbs, default_slurm
  
  # sets up default input.
  if outdir == None: outdir = getcwd() 

  which = "SNLCLUSTER" in environ
  if which: which = environ["SNLCLUSTER"] in ["redrock", "redmesa"]
  template = default_slurm if which else default_pbs

  # creates directory.
  directory = abspath( join(outdir, "pbs_scripts") )
  # saves pickle
  save(jobdict, join(outdir, "job_pickle"))
  # creates directory.
  with Changedir(directory) as pwd: pass 
  # gets runone 
  pyscript = __file__.replace(pathsplit(__file__)[1], "runone.py")
  # makes sure ppath is absolute.
  if ppath != None:
    kwargs["ppath"] = abspath(expanduser(ppath))
  # creates pbs script for each job.
  results = []
  for i, (job, name) in enumerate(jobdict.walk_through()):
    if job.is_tagged: continue
    mppwidth = mppalloc(job) if hasattr(mppalloc, "__call__") else mppalloc
    name = name.replace("/", ".")
    results.append( abspath(join(directory, name + ".pbs")) )
    with open(results[-1], "w") as file: 
      template( file, outdir=outdir, jobid=i, mppwidth=mppwidth, name=name,\
                pickle = join(outdir, "job_pickle"), pyscript=pyscript, **kwargs )
    print "wrote pbs script: %s." % (results[-1])
  return results





def fakerun(jobdict, outdir = None):
  """ Performs a fake run.

      Fake runs include *norun=True* as a job parameter. Whether this works or
      not depends on the functional. It is meant to create all directories and
      input file for a quick review of the input parameters.
      @param jobdict: otherwise fake runs the provided dictionary. 
      @type jobdict; L{JobDict} or None
      @param path: This will be the ouput directory.
  """
  from os import getcwd

  if outdir  == None: outdir = getcwd()
  for job, dirname in jobdict.walk_through(outdir):
    if not job.is_tagged: job.compute(outdir=dirname, norun=True)
