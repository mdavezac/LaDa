""" Submodule declaring the job-dictionary class. """
__docformat__ = "restructuredtext en"

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
      ``method(parameters)``), ``method`` is the callable object:

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

        >>> for name, job in jobtree.iteritems():
        >>>    result = job.compute(outdir=join(rootdir, name))


      The ``iteritems`` method does just that: it goes through the whole tree
      returning those branches where there is something to execute (eg
      ``job.functional != None``).
        
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
        >>> for name, job in jobtree.iteritems():
        >>>    job.compute(outdir=join(rootdir, name))

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
      ``jobdict`` saved to "pickle". This means that once all jobs are executed,
      ``bleed`` will find that ``jobdict`` is empty. To undo these changes and
      use the jobdict for, say, ouptut analysis, simply use ``iteritems``
      where possible (single pbs script or interactive
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
    try: self._functional = loads(string)
    except Exception as e:
      raise ValueError("Could not reload pickled functional. Caught Error:\n{0}".format(e))
  @functional.deleter
  def functional(self): self._functional = None

  @property
  def name(self):
     """ Returns the name of this dictionary, relative to root. """
     if self.parent == None: return "/"
     string = None
     for key, item in self.parent.children.iteritems():
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
    
        If the job does not exist, will create it.  A copy (copy.deepcopy) of
        value is inserted, rather than a simple shallow ref.
    """
    from re import split
    from copy import deepcopy
    from os.path import normpath, relpath, dirname, basename

    index = normpath(index)
    parentpath, childpath = dirname(index), basename(index)
    if len(parentpath) != 0: 
      assert parentpath in self, KeyError('Could not find parent job {0}.'.format(parentpath))
      mother = self[parentpath]
      parent = self.parent
      while parent != None:
        assert parent is not mother, KeyError('Will not set parent job of current job.')
    assert len(childpath) > 0 or childpath == '.', KeyError('Will not set current directory.')
    assert childpath != '..', KeyError('Will not set parent directory.')

    parent = self if len(parentpath) == 0 else self[parentpath]
    parent.children[childpath] = deepcopy(value)
    parent.children[childpath].parent = parent

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
    return sorted(self.children.iterkeys())
    
  def compute(self, **kwargs):
    """ Performs calculations over job list. """  

    if not self.is_job: return None
    kwargs.update(self.jobparams)
    return self.functional(**kwargs)

  def update(self, other, merge=False):
    """ Updates job and tree with other.
    
        :Parameters: 
           other : `JobDict`
             An other job tree from which to update.
           merge : bool
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
    for key, value in other.children.iteritems():
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
    for name in self.iterkeys():
      result += "  " + name + "\n"
    return result

  @property
  def untagged_jobs(self):
    """ Returns a string with only untagged jobs. """
    result = "Jobs: \n"
    for name, job in self.iteritems():
      if not job.is_tagged: result += "  " + name + "\n"
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
    from itertools import chain
    result = chain([u for u in self.__dict__ if u[0] != '_'], \
                   [u for u in dir(self.__class__) if u[0] != '_'], \
                   [u for u in self.jobparams.iterkeys() if u[0] != '_'])
    return list(set(result))

  def __getstate__(self):
    d = self.__dict__.copy()
    jobparams = d.pop("jobparams")
    return d, jobparams, "This is a JobDict pickle. Grep me!"
  def __setstate__(self, args):
    super(JobDict, self).__setattr__("jobparams", args[1])
    d = self.__dict__.update(args[0])

  def iteritems(self, outdir=''):
    """ Iterates over jobs. 

        :return: yields (directory, job):
          - directory is a suggested directory name with ``outdir`` as its root.
          - job is a Jobdict which contains something to execute.
    """
    from os.path import join
    # Yield this job if it exists.
    if self.is_job: yield outdir, self
    # Walk throught children jobdict.
    for name in self.subjobs():
      for u in self[name].iteritems(join(outdir, name)): yield u
  def itervalues(self): 
    """ Iterates over all jobs. """
    for name, job in self.iteritems(): yield job
  def iterkeys(self): 
    """ Iterates over all jobs. """
    for name, job in self.iteritems(): yield name
  def values(self):
    """ List of all jobs. """
    return [u for u in self.itervalues()]
  def keys(self):
    """ List of all jobs. """
    return [u for u in self.iterkeys()]
  def items(self):
    """ List of all jobs. """
    return [u for u in self.iteritems()]


  @property
  def nbjobs(self):
    """ Returns the number of jobs in tree. """
    return len([0 for j, o in self.iteritems() if not j.is_tagged])

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
    if index == '/': return True
    if index[0] == '/': return index[1:] in self.root
    names = split(r"(?<!\\)/", index) 
    if len(names) == 0: return False
    if len(names) == 1: return names[0] in self.children
    if names[0] not in self.children: return False
    new_index = normpath(index[len(names[0])+1:])
    if len(new_index) == 0: return True
    return new_index in self[names[0]]

  def __copy__(self):
    """ Performs a shallow copy of this job-dictionary.

        Shallow copies are made of all internal dictionaries children and
        jobparams. However, functional and jobparams values should the same
        object as self. The sub-branches of the returned dictionary are shallow
        copies of the sub-branches of self. In other words, the functional and
        refences in jobparams dictionary are in common between result and self,
        but nothing else.

        The returned dictionary does not have a parent!
    """
    from copy import copy
    # new job-dictionary.
    result = JobDict()
    result._functional = self._functional
    result.jobparams   = self.jobparams.copy()
    result.parent     = None
    for name, value in self.children.items():
      result.children[name] = copy(value)
      result.children[name].parent = result
    attrs = self.__dict__.copy()
    attrs.pop('jobparams')
    attrs.pop('parent')
    attrs.pop('children')
    attrs.pop('_functional')
    result.__dict__.update(attrs)
    return result


