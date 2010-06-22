""" Allows the creation of jobs. """

class JobDict(object):
  """ Holds all current jobs. """

  _pickle_filename = "job_pickle"
  """ Saves job to perform in current directory. """

  def __init__(self):
    super(JobDict, self).__init__()
    self.children = {}
    """ Jobs inheriting from this one. """
    self.job = {}
    """ This particular job. """
    self.parent = None
    """ Parent job. """
    
  def __getitem__(self, name): 
    """ Returns job description from the dictionary.

        If the value is a string, unpickles it before returning.
    """
    from re import split

    names = split("(?<!\\\)/", name)
    result = self
    for name in names[:-1]:
      if name not in result.children: 
        result.children[name] = JobDict()
        result.children[name].parent = result
      result = result.children[name]
    if names[-1] in result.children: 
      return result.children[names[-1]]
    elif names[-1] in result.job: return result.job[names[-1]]
    else:
      result.children[names[-1]] = JobDict()
      result.children[names[-1]].parent = result
      return result.children[names[-1]]

  def __setitem__(self, name, value): 
    """ Sets a job description in the dictionary.
    
        If value is not a string, then it is pickled and the
        result inserted into the job dictionary. Otherwise, it is expected the
        string is a pickle.
    """
    from copy import deepcopy
    from re import split

    names = split("(?<!\\)/", name)
    result = self
    for name in names[:-1]:
      assert name in result.children, RuntimeError("Job %s is not known." % (name))
      result = result.children[name]
    if names[-1] in result.children: 
      assert isinstance(value, JobDict),\
             RuntimeError("Cannot set job attribute %s. Job with same name exists." % (names[-1]))
      result.children[names[-1]]
    elif names[-1] in result.job:
      assert not isinstance(value, JobDict),\
             RuntimeError("Cannot set job %s. Job attribute with same name exists." % (names[-1]))
      result.children[names[-1]] = value
    if isinstance(value, JobDict): result.children[names[-1]] = deepcopy(value)
    else: result.job[names[-1]] = deepcopy(value)

  def has_job(self):
    """ True if self.job has keyword \"vasp\". """
    return "vasp" in self.job

  def items(self):
    """ Iterator over children jobs. """
    return self.children.items()
  def jobs(self):
    """ Iterator over children jobs. """
    return sorted(self.children.keys())
    
  def compute(self, **kwargs):
    """ Performs calculations over job list. """
    from copy import deepcopy

    kwargs.update(self.job)
    if "vasp" not in kwargs: return None
    vasp = kwargs.pop("vasp")
    args = kwargs.pop("args", ())
    assert hasattr(args, "__iter__"), RuntimeError("Functional argument \"args\" is not a sequence.")
    return vasp(*args, **kwargs)

  def update(self, other):
    """ Updates job and tree with other. """
    self.children.update(other.children)
    self.job.update(other.job)

current = JobDict()
""" Global with current joblist. """

def walk_through(jobdict = None, outdir = None):
  """ Generator to iterate over actual calculations. """
  from os import getcwd
  from os.path import join
  from copy import deepcopy

  if outdir == None: outdir = getcwd()
  if jobdict == None: jobdict = current
  # Yield this job if it exists.
  if jobdict.has_job():
    result = JobDict()
    result.job = jobdict.job
    yield result, outdir
  # Walk throught children jobdict.
  for name in jobdict.jobs():
    for u in walk_through(jobdict[name], join(outdir, name)): 
      yield u

def save(jobdict = None, path = None): 
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
  from os.path import exists
  from cPickle import load
  if path == None: path = "pickled_jobdict"
  if not exists(path): 
    print path, " does not exist."
    return
  print "Loading job list from", path, "."
  with open(path, "r") as file: return load(file)
