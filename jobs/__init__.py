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
__docformat__ = "restructuredtext en"
from ..opt.decorators import add_setter, broadcast_result, make_cached

__all__ = [ 'JobDict', 'walk_through', 'save', 'load', 'bleed', 'unbleed',\
            'unsucessfull', 'Extract' ]


class JobDict(object):
  """ Tree of jobs. 

      User Guide
      ==========

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


      The ``walkthrough`` method does just that: it goes through the whole tree
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

      ``bleed`` makes sure that once a job is accepted, it will not be accessed
      by any other process or pools of processes. ``bleed`` allows parallelization
      over pbs scripts and pools of MPI processes. However, it does change the
      C{jobdict} saved to "pickle". This means that once all jobs are executed,
      ``bleed`` (but not ``walk_through``) will find that C{jobdict} is empty. To
      undo these changes and use the jobdict for, say, ouptut analysis, simply
      use ``walk_through`` where possible (single pbs script or interactive
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
     



      Coding 
      ======

      JobDict has the following attributes:

        - children: A dict object holding instances of JobDict. These are the sub-jobs.
        - jobparams: All parameters regarding actual calculations. It contains,
          at start, only two predefined parameters.

           - functional: is the callable (preferably pickleable) to execute.
           - all others are keyword arguments.

        - functional: if None there are no jobs. Otherwise is a callable
          accepting jobparams as keyword arguments.
        - parent is an instance to the parent job (eg the instance which holds
          self in children) or None.
        - It may also have a _tagged attribute to check for bled/unbled jobs.

      The __getattr__, __setattr__, and __delattr__ have been rewired to
      perform on objects in jobparams. Note however that __setattr__ will not
      set new object in jobparams, but rather pass on the call to the parent
      class __setattr__. To set new job parameters, one should use
      ``add_param`` or jobparams directly.
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

  @property
  def functional(self):
    """ Returns current functional.
    
        The functional is implemented as a property to make sure that it is
        either None or a pickleable callable. Furthermore, a deepcopy of the functional is
        performed. This parameter can never be truly deleted.

        >>> del job.functional 

        is equivalent to:
        
        >>> job.functional = None
    """
    return self._functional

  @functional.setter
  def functional(self, value):
    from pickle import dumps, loads # ascertains pickle-ability, copies functional
    assert value == None or hasattr(value, "__call__"),\
           ValueError("job.functional should be either None(no job) or a callable.")
    # ascertains pickle-ability
    try: string = dumps(value)
    except Exception as e:
      raise ValueError("Could not pickle functional. Caught Error:\n{0}".format(e))
    else: self._functional = loads(string)
  @functional.deleter
  def functional(self): self._functional = None

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

  def update(self, other, merge=False):
    """ Updates job and tree with other.
    
        :Param other: An other job tree from which to update.
        :type other: `JobDict`
        :Param merge: 
          If false (default), then actual jobs in ``other`` completely
          overwrite actual jobs in ``self``. If False, then ``jobparams`` in
          ``self`` is updated with ``jobparams`` in ``other`` if either one is
          an actual job. If ``other`` is an actual job, then ``functional`` in
          ``self`` is overwritten. If ``other`` is not an actual job, then
          ``functional`` in ``self`` is not replaced.
         

        Updates the dictionaries of job parameters and sub-jobs. Actual jobs in
        ``other`` (eg with ``self.is_job=True``) will completely overwrite those in
        ``self``.  if items in ``other`` are found in ``self``, unless merge is
        set to true. This function is recurrent: subjobs are also updated.
    """
    for key, value in other.children.items():
      if key in self: self[key].update(value)
      else: self[key] = value

    if not merge:
      if not other.is_job: return
      self.jobparams = other.jobparams
      self.functional = other.functional
    else:
      if not (self.is_job or other.is_job): return
      self.jobparams.update(other.jobparams)
      if other.functional != None: self.functional = other.functional

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
    result = [u for u in self.__class__.__dict__ if u[0] != '_'] 
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

  def __contains__(self, index):
    """ Returns true if index a branch in the job-dictionary. """
    from re import split
    from os.path import normpath
    index = normpath(index)
    if index[0] == '/': return index[1:] in self.root
    names = split(r"(?<!\\)/", index) 
    if len(names) == 0: return False
    if len(names) == 1: return names[0] in self.children
    if names[0] not in self.children: return False
    new_index = normpath(index[len(names[0])+1:])
    if len(new_index) == 0: return True
    return new_index in self[names[0]]


def walk_through(jobdict, outdir = None, comm = None):
  """ Generator to iterate over actual calculations. 
      
      see ``JobDict`` description.
  """
  if isinstance(jobdict, str): jobdict = load(jobdict, comm=comm)
  for u in jobdict.walk_through(outdir): yield u

@broadcast_result(key=True)
def save(jobdict, path = None, overwrite=False, comm=None): 
  """ Pickles a job to file. 
 
      :keyword jobdict: A job-dictionary to pickle. 
      :type jobdict: `JobDict`
      :keyword path: 
          filename of file to which to save pickle. overwritten. If None then
          saves to "pickled_jobdict"
      :type path: str or None
      :keyword comm:
        Convenience parameter. Only root process actually saves.
        Other processes wait silently.
      :type comm: boost.mpi.comm
      :keyword overwrite: if True, then overwrites file.

      This method first acquire an exclusive lock on the file before writing
      (see `lada.opt.open_exclusive`).  This way not two processes can
      read/write to this file while using this function.
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
 
      :keyword path: Filename of a pickled jobdictionary.
      :keyword comm: MPI processes for which to read job-dictionary.
      :type comm: boost.mpi.communicator
      :return: Returns a JobDict object.

      This method first acquire an exclusive lock (using os dependent lockf) on
      the file before reading. This way not two processes can read/write to
      this file while using this function.
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

      :Parameters: 
        path 
          Filename of a pickled jobdictionary.
        outdir
          Root result directory. 
        comm
          Will broadcast yielded stuff from root. Because of file locking,
          this generator may freeze the system if not used correctly with mpi.

      :return: yields (job, directory), see ``walk_through``.
           - job: a job dictionary with the current job to execute.
           - directory: a suggested directory name with ``outdir`` as its root.

      This function alters the dictionary stored in `path`. If `path` is
      empty, then returns None, None. An exclusive lock is acquired before
      reading/writing to `path`.  This way, if using `bleed`, `save`,
      `load`, two processes will not step on each others jobs.
      This function is different from JobDict.bleed in that it alters a
      dictionary stored in a file.
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

      :Parameters:
        jobdict : JobDict or str
          a job dictionary instance or the path to a pickled job-dictionary.
        extractor
          Some instance with a directory attribute and a success attribute,
          capable of judging the success of an operation.

      :return: A JobDict instance with incomplete and unsuccessful jobs.

      This function will not modify the input dictionary.
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

class AbstractMassExtract(object): 
  """ Propagates extraction methods from different jobs. """

  def uncache(self): 
    """ Uncache values. """
    self.__dict__.pop("_cached_extractors", None)
    self.__dict__.pop("_cached_properties", None)


  def __init__(self, naked_end=True, unix_re=True, view=None, excludes=None):
    """ Initializes extraction object. 


        :Parameters:
          naked_end : bool
            True if should return value rather than dict when only one item.
          unix_re : bool
            converts regex patterns from unix-like expression.
          view : str or None
            Pattern which the job names must match to be included in the
            extraction.
          excludes : list of str or None
            List of patterns which the job names must *not* match to be
            included in the extraction.
    """
    from ..opt import RelativeDirectory

    super(AbstractMassExtract, self).__init__()

    self.naked_end = naked_end
    """ If True and dict to return contains only one item, returns value itself. """
    self.view = view
    """ The pattern which job-names should match. """
    self.unix_re = unix_re
    """ If True, then all regex matching is done using unix-command-line patterns. """
    self.excludes = excludes
    """ List of patterns to ignore. or None.

        ``self.unix_re`` determines whether these are unix-command-line like
        patterns or true python regex.
    """ 

  def _regex_pattern(self, pattern, flags=0):
    """ Returns a regular expression. """
    from ..opt import convert_from_unix_re
    return compile(pattern, flags) if not self.unix_re\
           else convert_from_unix_re(pattern)

  def walk_through(self):
    """ Generator to go through all relevant jobs. 
    
        :return: (name, extractor), where name is the name of the job, and
          extractor an extraction object.
    """
    abstract

  @property
  def view(self):
    """ A regex pattern which the name of extracted jobs should match.

        If None, then no match required. Should be a string, not an re object.
    """
    if self._view == None: return ""
    return self._view
  @view.setter
  def view(self, value): self._view = value

  def _extractors(self):
    """ Goes through all jobs and collects Extract if available. """
    from os.path import exists, join
    
    if "_cached_extractors" in self.__dict__:
      return self._cached_extractors
    result = {}
    for name, extract in self.walk_through(): result[name] = extract 
    self._cached_extractors = result
    return result

  def _regex_extractors(self):
    """ Loops through jobs in this view. """
    if self._view == "": 
      for key, value in self._extractors().items(): yield key, value
      return

    regex = self._regex_pattern(self.view)
    if self.excludes != None: excludes = [self._regex_pattern(u) for u in self.excludes]
    for key, value in self._extractors().items():
      if regex.match(key) == None: continue
      if self.excludes != None and any(u.match(key) != None for u in excludes): continue
      yield key, value

  def _properties(self): 
    """ Returns __dir__ special to the extraction itself. """
    if hasattr(self, "_cached_properties"): return self._cached_properties

    results = set([])
    for key, value in self._regex_extractors():
      results |= set([u for u in dir(value) if u[0] != '_'])
    self._cached_properties = results
    return results

  def __dir__(self): 
    results =   set([u for u in self.__dict__ if u[0] != '_']) \
              | set([u for u in self.__class__.__dict__ if u[0] != '_']) \
              | self._properties()
    return list(results)

  def __getattr__(self, name): 
    """ Returns extracted values. """
    assert name not in ["_cached_extractors", "_cached_properties"],\
           AttributeError("Unknown attribute {0}.".format(name))
    assert name in self._properties(), AttributeError("Unknown attribute {0}.".format(name))

    result = {}
    for key, value in self._regex_extractors():
      try: result[key] = getattr(value, name)
      except: result.pop(key, None)
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return result

  def __getitem__(self, name):
    """ Returns a view of the current job-dictionary. """
    from os.path import normpath, join
    if name[0] == '/': return self.copy(view=name)
    path = normpath(join('/', join(self.view, name)))
    return self.copy(view=path)

  @property
  def jobs(self):
    """ List of jobnames (satisfying the current view). """
    return [key for key, value in self._regex_extractors()]

  @property
  def children(self):
    """ next set of minimal regex. """
    from os.path import join, normpath
    regex = self._regex_pattern(self.view)

    jobs = self.jobs
    if len(self.jobs) < 2: return 
    children = set()
    if len(self.view) == 0:
      for name in self.jobs:
        children.add(name[:1+name[1:].find('/')])
    else:
      for name in self.jobs:
        where = regex.match(name)
        if len(name) == where.end() +1: continue
        first_index = name[where.end():].find('/')
        if first_index == -1: continue
        first_index += where.end() + 1
        if first_index >= len(name): continue
        end = name[first_index:].find('/')
        if end == -1: children.add(name)
        else: children.add(name[:end + first_index])
    
    for child in children: yield self.copy(view=child)

  def grep(self, regex, flags=0, yield_match=False):
    """ Yields views for children with fullnames matching the regex.
    
        :Parameters:
          regex : str
            The regular expression which the fullnames should match. Whether
            this is a python regex, or something which behaves like the unix
            command-line depends on ``self.unix_re``.
          flags : int
             Flags from ``re`` to use when compilling the regex pattern.
          yield match: bool
             If True, will yield a two tuple, where the second item is the
             match object.
             If False, only the view is yielded.
             This option is not available (or meaningfull) if ``self.unix_re``
             is True.

        The match is successful if the regex is matched using python's
        `re.search <http://docs.python.org/library/re.html#re.search>`_ method.

        Only the innermost view of each match is given. In other words, if a
        view is yielded, its subviews will not be yielded.

        If the current view matches the regex, then it alone is yielded. 
    """
    assert not (yield_match and self.unix_re),\
           ValueError("unix_re and yield_matc cannot be both true.") 
    reg = self._regex_pattern(regex, flags)

    found = reg.search(self.view)
    if found != None and yield_match:       yield self; return
    elif found != None and not yield_match: yield self, found; return
    
    for child in self.children:
      found = reg.search(self.view)
      if reg.search(child.view) == None:# goes to next level. 
        for grandchild in child.grep(regex, flags, yield_match): yield grandchild
      elif yield_match: yield child, found
      else: yield child

  def __getstate__(self):
    d = self.__dict__.copy()
    d.pop("comm", None)
    if "_rootdir" in d: d["_rootdir"].hook = None
    return d

  def __setstate__(self, arg):
    self.__dict__.update(arg)
    self.comm = None
    if "_rootdir" in d: d["_rootdir"].hook = self.uncache
       
  def solo(self):
    """ Extraction on a single process.
  
        Sometimes, it is practical to perform extractions on a single process
        only, eg without blocking mpi calls. C{self.``solo}()`` returns an
        extractor for a single process:
        
        >>> # prints only on proc 0.
        >>> if boost.mpi.world.rank == 0: print extract.solo().structure
    """
    if self.comm == None: return self

    from copy import deepcopy
    copy = deepcopy(self)
    return copy

  def __copy__(self):
    """ Returns a shallow copy. """
    result = self.__class__(self.rootdir)
    result.__dict__.update(self.__dict__)
    return result

  def copy(self, **kwargs):
    """ Returns a shallow copy. 
    
        :Param kwargs:  Any keyword attribute will modify the corresponding
          attribute of the copy.
    """
    from copy import copy
    result = copy(self)
    for key, value in kwargs.items(): setattr(result, key, value)
    return result



class MassExtract(AbstractMassExtract): 
  """ Propagates extraction methods from different jobs. 
  
      Collects extractors across all jobs (for which job.functional.Extract
      exist). The results are presented as attributes of an instance of
      MassExtract, and arranged as directory where the key is the name of the
      job and the value obtained from an instance of that job's Extract. This
      class is set-up to fail silently, and hence is of limited use for
      diagnosis.
  """

  def __init__(self, path=None, naked_end=True, unix_re=True,\
                     view=None, excludes=None, comm=None):
    """ Initializes extraction object. 
 
        :Parameters:
          path : str or None
            Pickled jobdictioanary for which to extract stuff. If None, will
            attempt to use the current jobdictionary.
          naked_end : bool
            True if should return value rather than dict when only one item.
          unix_re : bool
            converts regex patterns from unix-like expression.
          view : str or None
            Pattern which the job names must match to be included in the
            extraction.
          excludes : list of str or None
            List of patterns which the job names must *not* match to be
            included in the extraction.
          comm : boost.mpi.communicator
             All processes will be syncronized.
    """
    super(MassExtract, self).__init__( naked_end=naked_end, unix_re=unix_re, \
                                       view=view, excludes=excludes)

    from os.path import isdir, isfile, exists, dirname, abspath

    self.rootdir = path # Path to the job dictionary.
    self.comm = comm
    """ Communicator if any. """

    if path != None: self._extractors() # gets stuff cached.

  @property
  def view(self):
    """ A regex pattern which the name of extracted jobs should match.

        If None, then no match required. Should be a string, not an re object.
    """
    if self._view == None:
      try: from IPython.ipapi import get as get_ipy
      except ImportError: raise AttributeError("path not set.")
      ip = get_ipy()
      if "current_jobdict" not in ip.user_ns:
        print "No current jobdictionary."
        return
      return ip.user_ns["current_jobdict"].name
      return 
    return self._view
  @view.setter
  def view(self, value): self._view = value

  @property
  def rootdir(self): 
    """ Root directory of the jobdictionary. """
    from os.path import dirname

    if self._rootdir == None: 
      try: from IPython.ipapi import get as get_ipy
      except ImportError: raise AttributeError("path not set.")
      ip = get_ipy()
      if "current_jobdict_path" not in ip.user_ns:
        print "No current jobdictionary path."
        return
      return dirname(ip.user_ns["current_jobdict_path"])

    return dirname(self._rootdir.path)
  @rootdir.setter
  def rootdir(self, value):
    from ..opt import RelativeDirectory
    if value == None:
      self._rootdir = None
      return
    self._rootdir = RelativeDirectory(value, hook=self.uncache)
    del self._jobdict
  @rootdir.deleter
  def rootdir(self): self._rootdir = None

  @property
  def jobdict(self):
    if self._rootdir == None: 
      try: from IPython.ipapi import get as get_ipy
      except ImportError: raise AttributeError("path not set.")
      ip = get_ipy()
      if "current_jobdict" not in ip.user_ns:
        print "No current jobdictionary."
        return
      return ip.user_ns["current_jobdict"].root
    if "_jobdict" not in self.__dict__: self._jobdict = load(self._rootdir.path, self.comm)
    return self._jobdict.root

  def walk_through(self):
    """ Generator to go through all relevant jobs.  
    
        :return: (name, extractor), where name is the name of the job, and
          extractor an extraction object.
    """
    from os.path import exists, join
    
    for job, name in self.jobdict.walk_through():
      if job.is_tagged: continue
      if not hasattr(job.functional, "Extract"): continue
      try: extract = job.functional.Extract(join(self.rootdir, name), comm = self.comm)
      except: pass
      else: yield job.name, extract


class JobParams(AbstractMassExtract):
  """ Get and sets job parameters for a job-dictionary. """
  def __init__( self, jobdict = None, only_existing=True,  naked_end=True,\
                unix_re=True, view = None, excludes = None):
    """ Initializes job-parameters.

        :Parameters:
          jobdict : None or JobDict
	    The jobdictionary for which to get/set parameters. If None, will
            look for ipython's current_jobdict.
          only_existing : bool
	    If True (default), will never create new parameters. It will 
            modify only existing parameters.
          naked_end : bool
	    If True, if the returned dictionary contains only one item, and if
	    that item corresponds to the root of the jobdictionary being
            explored, returns that item alone.
          unix_re : bool
            converts regex patterns from unix-like expression.
          view : None or str
            Grepable to examine.
    """
    super(JobParams, self).__init__( naked_end=naked_end, unix_re=unix_re, \
                                     view=view, excludes=excludes)

    super(JobParams, self).__setattr__("_jobdict", None)
    self._jobdict = jobdict
    """ Job-dictionary for which to get/set parameters. """
    super(JobParams, self).__setattr__("only_existing", None)
    self.only_existing = only_existing
    """ Only modifies parameter which already exist. """

  @property
  def jobdict(self):
    """ Jobdictionary for which to get/set parameters. """
    if self._jobdict == None:
      try: from IPython.ipapi import get as get_ipy
      except ImportError: raise AttributeError("jobdict not set.")
      else:
        ip = get_ipy()
        if "current_jobdict" not in ip.user_ns:
          print "No current jobdictionary."
          return
        return ip.user_ns["current_jobdict"].root
    return self._jobdict.root
  @jobdict.setter
  def jobdict(self, value): self._jobdict = value

  @property
  def view(self):
    """ A regex pattern which the name of extracted jobs should match.

        If None, then no match required. Should be a string, not an re object.
    """
    if self._view == None:
      try: from IPython.ipapi import get as get_ipy
      except ImportError: raise AttributeError("path not set.")
      ip = get_ipy()
      if "current_jobdict" not in ip.user_ns:
        print "No current jobdictionary."
        return
      return ip.user_ns["current_jobdict"].name
      return 
    return self._view
  @view.setter
  def view(self, value): self._view = value

  def walk_through(self):
    """ Loops through all correct jobs. """
    for job, name in self.jobdict.walk_through():
      if not job.is_tagged: yield job.name, job

  def __setattr__(self, name, value):
    """ Returns dictionary with job parameters for each job. """
    from re import match
    # initialization not done yet.
    if "only_existing" not in self.__dict__: super(JobParams, self).__setattr__(name, value)
    # some cached attributes.
    if name == "_cached_extractors" or name == "_cached_properties":
      super(JobParams, self).__setattr__(name, value)
    # other cached attributes.
    if match("_cached_attr\S+", name): super(JobParams, self).__setattr__(name, value)
    try: super(JobParams, self).__getattribute__(name)
    except AttributeError: 
      for jobname, job in self._regex_extractors():
        if hasattr(job, name): setattr(job, name, value)
        elif not self.only_existing: job.jobparams[name] = value
    else: super(JobParams, self).__setattr__(name, value)

  def __delattr__(self, name):
    try: super(JobParams, self).__getattribute__(name)
    except AttributeError: 
      for jobname, job in self._regex_extractors():
        if hasattr(job, name): delattr(job, name)
    else: super(JobParams, self).__delattr__(name, value)

  def _properties(self):
    """ Attributes which already exist. """
    result = set()
    for name, job in self._regex_extractors(): result |= set(job.jobparams.keys())
    return result
 
