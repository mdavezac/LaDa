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
from ..opt import RelativeDirectory
from ..opt.decorators import add_setter, broadcast_result, make_cached

__all__ = [ 'JobDict', 'walk_through', 'save', 'load', 'bleed', 'unbleed',\
            'unsucessfull', 'Extract' ]


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
      "python pickle"), then method should be pickleable. All arguments to the
      functional are specified as keyword arguments: 

        >>> subjob1.jobparams["keywordname"] = value

      Once they have been specified, keyword arguments can be accessed directly:

        >>> # first specify
        >>> subjob1.jobparams["keywordname"] = value1
        >>> # then change
        >>> subjob1.keywordname = value2

      Where keywordname in the last line should be substituted for its actual
      value (eg structure, nelect...). Note that there is no options to insert
      non-keyword argument. This is a feature and most likely will stay that
      way.

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
        >>> subjob1.add_param = "whatever", value1
        >>> subjob1.whatever = value2
        >>> # actually performs calculations.
        >>> for job, name in walkthrough(jobdict=jobtree, outdir="root_result_directory"):
        >>>    job.compute(outdir=name)

      Using the jobtree defined above, this would end-up calling:

        >>> metod(whatever=value1, outdir="root_result_directory/subjoname1") 


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
           - all others are keyword arguments.
        - parent is an instance to the parent job (eg the instance which holds
          self in children) or None.
        - It may also have a _tagged attribute to check for bled/unbled jobs.
      The __getattr__, __setattr__, and __delattr__ have been rewired to
      perform on objects in jobparams. Note however that __setattr__ will not
      set new object in jobparams, but rather pass on the call to the parent
      class __setattr__. To set new job parameters, one should use
      L{add_param} or jobparams directly.
  """

  def __init__(self):
    super(JobDict, self).__init__()
    # List of subjobs (as in subdirectories). 
    super(JobDict, self).__setattr__("children", {})
    # This particular job. 
    super(JobDict, self).__setattr__("jobparams", {})
    # This particular job is not set. 
    super(JobDict, self).__setattr__("_functional", None)
    # Parent job. 
    super(JobDict, self).__setattr__("parent", None)

  def _get_functional(self):
    """ Returns current functional.
    
        The functional is implemented as a property to make sure that it is
        either None or a pickleable callable. Furthermore, a deepcopy of the functional is
        performed. This parameter can never be truly deleted.
          >>> del job.functional 
        is equivalent to:
          >>> job.functional = None
    """
    return self._functional
  def _set_functional(self, value):
    from pickle import dumps, loads # ascertains pickle-ability, copies functional
    assert value == None or hasattr(value, "__call__"),\
           ValueError("job.functional should be either None(no job) or a callable.")
    # ascertains pickle-ability
    try: string = dumps(value)
    except Exception as e:
      raise ValueError("Could not pickle functional. Caught Error:\n{0}".format(e))
    else: self._functional = loads(string)
  def _del_functional(self): self._functional = None

  functional = property(_get_functional, _set_functional, _del_functional)
    
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
      else: raise KeyError("job " + index + " does not exist.")
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
    raise KeyError("job " + index + " does not exist.")

  def __setitem__(self, index, value): 
    """ Sets job/subjob description in the dictionary.
    
        If the job does not exist, will create it.
        A copy (copy.deepcopy) of value is inserted, rather than a simple
        shallow ref.
    """
    from re import split
    from copy import deepcopy
    from os.path import normpath, relpath

    index = normpath(index)
    assert index not in ["", ".", None], KeyError("Will not set self.")
    assert index[0], KeyError("Will not set root: " + index + ".")
    assert isinstance(value, JobDict), \
           ValueError("Only JobDict instances can be used a job dictionaries.")

    result = self.__div__(index+"/..")
    name = relpath(index, index+"/..")
    result.children[name] = deepcopy(value)
    result.children[name].parent = result

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
    """ True if functional is not None. """
    return self.functional != None

  def subjobs(self):
    """ Iterator over children jobs. """
    return sorted(self.children.keys())
    
  def compute(self, **kwargs):
    """ Performs calculations over job list. """  

    if not self.is_job: return None
    kwargs.update(self.jobparams)
    return self.functional(**kwargs)

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
    if name in self.__dict__: return self.__dict__.pop(name)
    if name in self.jobparams: return self.jobparams.pop(name)
    raise AttributeError("Unknown job attribute " + name + ".")

  def __getattr__(self, name):
    """ Returns job attribute. """
    if name in self.jobparams: return self.jobparams[name]
    raise AttributeError("Unknown job attribute " + name + ".")

  def __setattr__(self, name, value):
    """ Sets job attribute. """
    from pickle import dumps
    if name in self.jobparams:
      try: dumps(value)
      except Exception as e:
        raise ValueError("Could not pickle job-parameter. Caught error:\n{0}".format(e))
      else: self.jobparams[name] = value
    else: super(JobDict, self).__setattr__(name, value)

  def __dir__(self):
    result = [u for u in self.__dict__ if u[0] != '_'] 
    result.extend([u for u in self.jobparams.keys() if u[0] != '_'])
    return list(set(result))

  def __getstate__(self):
    d = self.__dict__.copy()
    jobparams = d.pop("jobparams")
    return d, jobparams, "This is a JobDict pickle. Grep me!"
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


class Extract(object): 
  """ Propagates extraction methods from different jobs. """
  root = RelativeDirectory()
  """ Root directory of the job. """

  def __init__(self, path, jobdict=None, comm = None):
    """ Initializes extraction object. 


        :Parameters:
          - `path` : If `jobdict` is None, then should point to a pickled
             dictionary. If `jobdict` is a JobDict instance, then it should
             point the directory where calculations are saved.
          - `jobdict` : None, or a jobdictionary. If it is None, then `path`
             should point to a pickled job-dictionary. If it is a jobdictonary,
             then `path` should point to a directory where calculations where performed.
          - `comm` : an boost.mpi.communicator instance.
    """
    from os.path import isdir, isfile, exists, dirname
    assert exists(path), IOErrror("{0} does not exist.".format(path))
    if jobdict == None:
      assert isfile(path), IOError("{0} is not a file.".format(path))
      with open(path, "r") as file: self.jobdict = load(path, comm)
      self.root = dirname(path)
    else: 
      assert isdir(path), IOError("{0} is not a directory.".format(path))
      self.root = path
      self.jobdict = jobdict
    self.comm = comm

  def _extractors(self):
    """ Goes through all jobs and collects Extract if available. """
    from os.path import exists, join
    
    if hasattr(self, "_cached_extractors"): return self._cached_extractors
    result = {}
    for job, name in self.jobdict.walk_through():
      if not hasattr(job.is_tagged, "Extract"): continue
      if job.is_tagged: continue
      if not exists(join(self.root, name)): continue
      try: result[name] = job.functional.Extract(join(self.root, name), comm = self.comm)
      except: pass
    self._cached_extractors = result
    return result

  def _properties(self): 
    """ Returns cached __dir__ result. """
    if hasattr(self, "_cached_properties"): return self._cached_properties
    results = set([])
    for key, value in self._extractors().items():
      for name in dir(value):
        if name[0] == '_': continue
        results.add(name)
    self._cached_properties = results
    return results


  def __dir__(self): 
    results = set([u for u in self.__dict__ if u[0] != '_']) | self._properties()
    return list(results)

  def __getattr__(self, name): 
    """ Returns extracted values. """
    if name == "_cached_extractors" or name == "_cached_properties": 
      raise AttributeError("Unknown attribute {0}.".format(name))
    if name in self._properties(): 
      result = {}
      for key, value in self._extractors().items():
        try: result[key] = getattr(value, name)
        except: result.pop(key, None)
      return result
    raise AttributeError("Unknown attribute {0}.".format(name))


  def uncache(self): 
    """ Uncache values. """
    self.__dict__.pop("_cached_extractors", None)
    self.__dict__.pop("_cached_properties", None)

  def __getstate__(self):
    d = self.__dict__.copy()
    d.pop("comm", None)
    return d

  def __setstate__(self, arg):
    self.__dict__.update(arg)
    self.comm = None
       
     



    
