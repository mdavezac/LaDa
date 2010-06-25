""" Allows the creation of jobs. 

    Contains the following methods and classes.
       - JobDict: tree of calculations and subcalculations.
       - walkthrough: iterates over calculation in a tree.
       - save: pickles a tree to file.
       - load: loads a pickled tree from file.
       - pbs_script: Creates pbs-script to perform calculations on a tree.
           Calculations can be parallelized simultaneously over different pbs
           scripts and over different pools of processes within each pbs script.
"""
from ..opt.decorators import add_setter

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
    assert index != "." 
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

        see L{JobDict} description.
    """
    for u in walk_through(jobdict=self, outdir=outdir): yield u

current = JobDict()
""" Global with current joblist. """

def walk_through(jobdict = None, outdir = None):
  """ Generator to iterate over actual calculations. 
      
      see L{JobDict} description.
  """
  from os.path import join

  if outdir == None: outdir = ""
  if jobdict == None: jobdict = current
  # Yield this job if it exists.
  if jobdict.has_job(): yield jobdict, outdir
  # Walk throught children jobdict.
  for name in jobdict.subjobs():
    for u in walk_through(jobdict[name], join(outdir, name)): 
      yield u

def save(jobdict = None, path = None): 
  """ Pickles a job to file.

      @param jobdict: A jobtree to pickle. If None, uses L{jobs.current}.
      @type jobdict: JobDict
      @param path: filename of file to which to save pickle. overwritten. If
        None then saves to "pickled_jobdict"
  """ 
  from os.path import exists
  from cPickle import dump
  if path == None: path = "pickled_jobdict"
  if jobdict == None: jobdict = current
  if exists(path): 
    print path, "exists. Please delete first if you want to save the job dictionary."
    return
  with open(path, "w") as file: dump(jobdict, file)
  print "Saved job dictionary to %s." % (path)

def load(path = None): 
  """ Unpickles a job from file.

      @param path: filename from which to load pickle. 
        If None then saves to "pickled_jobdict"
      @return: Returns a JobDict object.
  """ 
  from os.path import exists
  from cPickle import load
  if path == None: path = "pickled_jobdict"
  assert exists(path), IOError("File " + path + " does not exist.")
  print "Loading job list from", path, "."
  with open(path, "r") as file: return load(file)


def pbs_scripts( outdir = None, jobdict = None, template = None, \
                 pbspools = 1, pickle = None, **kwargs):
  """ Parallelizes jobs over different pbs scrits. 

      The root directory (outdir) is created if it does not exist, and the
      pbs-scripts placed there. The pbs-scripts are not launched. 

      @param outdir: root directory of calculation. Current working directory if None.
      @param jobdict: job-tree instance over which to parallelize.
      @type jobdict: L{JobDict}
      @param template: PBS-script template to use. See L{jobs.templates}. Default: L{jobs.templates.default_pbs}.
      @type template: callable.
      @param pbspools: Number of pbs scripts to issue. 
      @param pickle: Name of file too which L{JobDict} will be pickled.
      @param kwargs: Passed on to template.
      @return: List of pbs-script filenames.
  """
  from os import getcwd, makedirs
  from os.path import abspath, join, exists
  from ..opt.changedir import Changedir
  from templates import default_pbs 
 
  if jobdict == None: jobdict = current
  if template == None: template = default_pbs
  outdir = getcwd() if outdir == None else abspath(outdir)
  if pickle == None: pickle = "job_pickle"

  # creates result dictionary.
  if not exists(outdir): makedirs(outdir)
  # pickle jobs.
  save(jobdict = jobdict, path = join(outdir, pickle))
  if pbspools <= 1:
    filpath = join(outdir, "automatic_pbs_job")
    with open(filepath, "w") as file:
      template(file, pbspools=pbspools, pickle=join(outdir, pickle), **kwargs)
    return [filepath]
  else:
    result = []
    for p in  range(pbspools+1):
      filepath = join(outdir, "automatic_pbs_job_%i" % (p))
      with open(filepath, "w") as file:
        template(file, pbspools=pbspools, pickle=join(outdir, pickle), **kwargs)
      result.append(filepath)

  return result

