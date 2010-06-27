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
           Calculations can be parallelized simultaneously over different pbs
           scripts and over different pools of processes within each pbs script.
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

        >>> subjob1.vasp = some callable object.

      Where a callable object is anything which can be called (eg, as in
      C{method(parameters)}), C{method} is the callable object:

        >>> subjob1.vasp = method

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

      Having a jobtree is great executing it is better:

        >>> for job, name in jobtree.walkthrough("root_result_directory"):
        >>>    job.compute(outdir=name)

      The L{walkthrough} method does just that: it goes through the whole tree
      returning those branches where there is something to execute (eg
      C{job.vasp != None}). It also returns outdir which is the subdirectory
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
        >>> subjob1.vasp = method
        >>> subjob1.args = (arg1, arg2) # two arguments.
        >>> subjob1.add_param = "whatever", value1
        >>> subjob1.whatever = value2
        >>> # actually performs calculations.
        >>> for job, name in walkthrough(jobdict=jobtree, outdir="root_result_directory"):
        >>>    job.compute(outdir=name)

      Using the jobtree defined above, this would end-up calling:

        >>> metod(arg1, arg2, whatever=value1, outdir="root_result_directory/subjoname1") 

      Coding: JobDict has name attributes:
        - children: A dict object holding instances of JobDict. These are the sub-jobs.
        - jobparams: All parameters regarding actual calculations. It contains,
          at start, only two predefined parameters.
           - vasp: is the callable (preferably pickleable) to execute.
           - args: is a tuple of arguments (eg. C{vasp(*args)}
           - all others are keyword arguments.
        - parent is an instance to the parent job (eg the instance which holds
          self in children) or None.
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
    self.jobparams["vasp"] = None
    # no arguments yet.
    self.jobparams["args"] = (None)
    
  def __getitem__(self, index): 
    """ Returns job description from the dictionary.

        If the job does not exist, will create it.
    """
    from re import split
    from os.path import normpath

    index = normpath(index)
    if index == "" or index == None or index == ".": return self
    if index[0] == "/":  # could create infinit loop.
      result = self
      while result.parent != None: result = result.parent
      return result[index[1:]]

    result = self
    names = split(r"(?<!\\)/", index)
    for i, name in enumerate(names):
      if name == "..":
        assert result.parent != None, RuntimeError("Cannot go below root level.")
        result = result.parent
      if name in result.children: result = result.children[name]
      elif name in result.jobparams:
        assert i+1 == len(names), KeyError("job or job parameter " + index + " does not exist.") 
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
      assert id(self) != id(deletee), RuntimeError("Will not commit suicide.")
      parent = self.parent
      while parent != None: 
        assert id(parent) != id(deletee), RuntimeError("Will not go Oedipus on you.")
        parent = parent.parent

    parent = self[index+"/.."]
    name = relpath(index, index+"/..")
    if name in parent.children:
      assert id(self) != id(parent.children[name]),\
             RuntimeError("Will not delete self.")
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
    assert index not in ["", ".", None], RuntimeError("Will not set self.")
    assert index[0], RuntimeError("Will not set root: " + index + ".")

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

  def has_job(self):
    """ True if self.job has keyword \"vasp\". """
    if "vasp" not in self.jobparams: return False
    return self.jobparams["vasp"] != None

  def items(self):
    """ Iterator over children jobs. """
    return self.children.items()
  def subjobs(self):
    """ Iterator over children jobs. """
    return sorted(self.children.keys())
    
  def compute(self, **kwargs):
    """ Performs calculations over job list. """

    kwargs.update(self.jobparams)
    if "vasp" not in kwargs: return None
    vasp = kwargs.pop("vasp")
    if vasp == None: return
    args = kwargs.pop("args", ())
    if args == None: args = ()
    assert hasattr(args, "__iter__"), RuntimeError("Functional argument \"args\" is not a sequence.")
    return vasp(*args, **kwargs)

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

 
  def __delattr__(self, name):
    """ Deletes job attribute. """
    if name in self.jobparams: return self.jobparams.pop(name)
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
        The only preset job parameter are "vasp" and "args", as stand in for
        the functional to use and the tuple to use as its arguments.
    """
    assert len(args) < 3, ValueError("Too many parameters given to add_param: %s." % (args))
    name = args[0]
    value = None if len(args) == 1 else args[1]
    self.jobparams[name] = value

  def __getstate__(self):
    d = self.__dict__.copy()
    params = d.pop("jobparams")
    return d, params
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
    if self.has_job(): yield self, outdir
    # Walk throught children jobdict.
    for name in self.subjobs():
      for u in self[name].walk_through(join(outdir, name)): 
        yield u

  @property
  def nbjobs(self):
    """ Returns the number of jobs in tree. """
    return len([u for u in self.walk_through()])

  def pop_first(self):
    """ Retrieves and removes first actual job. 
    
        The first jobs jobparam is reset. Subjobs are not, of course.
        Returns None if no jobs.
    """
    for job, outdir in self.walk_through():
      # new dictionary
      result = JobDict()
      # swap functional stuff
      result.jobparams, job.jobparams = job.jobparams, result.jobparams
      # adds subjobs.
      result.children = job.children
      return result, outdir
    return None
      

current = JobDict()
""" Global with current joblist. """

def walk_through(jobdict = None, outdir = None):
  """ Generator to iterate over actual calculations. 
      
      see L{JobDict} description.
  """
  if jobdict == None: jobdict = current
  for u in jobdict.walk_through("outdir"): yield u

@broadcast_result(key=True)
def save(jobdict = None, path = None, overwrite=False, comm=None): 
  """ Pickles a job to file.

      This method first acquire an exclusive lock (using os dependent flock) on
      the file before writing. This way not two processes can read/write to
      this file while using this function.
      @param jobdict: A jobtree to pickle. If None, uses L{jobs.current}.
      @type jobdict: JobDict
      @param path: filename of file to which to save pickle. overwritten. If
        None then saves to "pickled_jobdict"
      @param comm: Only root process gets to do anything.
      @type comm: boost.mpi.communicator
      @param overwrite: if True, then overwrites file.
  """ 
  from os.path import exists
  from cPickle import dump
  from ..opt import open_exclusive
  if path == None: path = "pickled_jobdict"
  if jobdict == None: jobdict = current
  if exists(path) and not overwrite: 
    print path, "exists. Please delete first if you want to save the job dictionary."
    return
  with open_exclusive(path, "w") as file: dump(jobdict, file)
  print "Saved job dictionary to %s." % (path)

@broadcast_result(key=True)
def load(path = None, comm = None): 
  """ Unpickles a job from file.

      This method first acquire an exclusive lock (using os dependent flock) on
      the file before reading. This way not two processes can read/write to
      this file while using this function.
      @param path: filename from which to load pickle. 
        If None then saves to "pickled_jobdict"
      @param comm: Broadcasts from root process.
      @type comm: boost.mpi.communicator
      @return: Returns a JobDict object.
  """ 
  from fcntl import flock, LOCK_EX, LOCK_UN
  from os.path import exists
  from cPickle import load as load_pickle
  from ..opt import open_exclusive
  if path == None: path = "pickled_jobdict"
  assert exists(path), IOError("File " + path + " does not exist.")
  print "Loading job list from", path, "."
  with open_exclusive(path, "r") as file: return load_pickle(file)

def bleed(path=None, outdir=None, comm=None): 
  """ Generator which deepletes a job dictionary of its jobs. 

      This function alters the dictionary path. If C{path} is empty, then
      returns None, None. An exclusive lock is acquired before reading/writing
      to C{path}.  This way, if using L{bleed}, L{save}, L{load}, two processes
      will not step on each others jobs.
      @param: Filename of a pickled jobdictionary.
      @outdir: Root result directory. 
      @comm: Will broadcast yielded stuff from root. Because of file locking,
             this generator may freeze the system if not used correctly with mpi.
      @return: yields (job, directory), see L{walk_through}.
           - job: a job dictionary with the current job to execute.
           - directory: a suggested directory name with L{outdir} as its root.
  """
  from fcntl import flock, LOCK_EX, LOCK_UN
  from os.path import join, exists
  from cPickle import load as load_pickle, dump
  from ..opt import open_exclusive
  from boost.mpi import broadcast
  if path == None: path = "pickled_jobdict"

  is_root = True if comm == None else comm.rank == 0
  if is_root:
    while True:
      if not exists(path): 
        print "Job dictionary", path, "does not exist."
        return
      with open_exclusive(path, "r") as file:
        # tries to load file
        try: jobdict = load_pickle(file)
        except EOFError: break
        # Checks if there are any jobs.
        if jobdict.nbjobs == 0: break
        # Pops first job.
        job, directory = jobdict.pop_first()
        # writes modified dictionary to path.
        with open(path, "w") as newfile: dump(jobdict, newfile)
      broadcast(comm, (job, join(outdir, directory)), 0)
      yield job, join(outdir, directory)
    broadcast(comm, None, 0)
    return
  else: 
    while True:
      result =  broadcast(comm, root=0)
      if result == None: return
      else: yield result
  


def pbs_scripts( outdir = None, jobdict = None, template = None, pbspools = 1,\
                 name = None, pickle = None, pyscript = None, **kwargs):
  """ Parallelizes jobs over different pbs scrits. 

      The root directory (outdir) is created if it does not exist, and the
      pbs-scripts placed there. The pbs-scripts are not launched. 

      @param outdir: root directory of calculation. Current working directory if None.
      @param jobdict: job-tree instance over which to parallelize.
      @type jobdict: L{JobDict}
      @param template: PBS-script template to use. See L{jobs.templates}.
        Default: L{jobs.templates.default_pbs}.
      @type template: callable.
      @param pbspools: Number of pbs scripts to issue. If 0 or negative, will
        issue one script per job.
      @param pickle: Filename to which L{JobDict} will be pickled.
      @param pyscript: Filename of python script to execute. Copied to outdir.
      @param kwargs: Passed on to template.
      @return: List of pbs-script filenames.
  """
  from shutil import copy
  from os import getcwd, makedirs
  from os.path import abspath, join, exists, relpath, samefile, split as pathsplit
  from ..opt.changedir import Changedir
  from templates import default_pbs 
 
  if jobdict == None: jobdict = current
  if template == None: template = default_pbs
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


  if pbspools > jobdict.nbjobs or pbspools < 1: pbspools = jobdict.nbjobs
  if "procpools" in kwargs:
    assert kwargs["procpools"] * pbspools <= jobdict.nbjobs, \
        ValueError( "Requested for more ressources than there are jobs:\n"\
                    "  - Number of jobs: %i\n"\
                    "  - Number of pbs scripts: %i\n"\
                    "  - Number of process pools: %i\n"\
                    % (jobdict.nbjobs, pbspools, kwargs["procpools"]) )

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
  for p in  range(pbspools+1):
    # writes pbs script.
    filepath = join(outdir, "launchme")
    jobname = name
    if pbspools > 1: 
      filepath += "_" + str(p)
      jobname += "-" + str(p)
    with open(filepath+".pbs", "w") as file:
      template( file, pbspools=pbspools, npbs=p, pickle=pickle, outdir=abspath(outdir),\
                name=jobname, pyscript=pyscript_filename, **kwargs)
    result.append(filepath)

  return result

