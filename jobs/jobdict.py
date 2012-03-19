""" Submodule declaring the job-dictionary class. """
__docformat__ = "restructuredtext en"

class JobDict(object):
  """ Tree/dictionary of jobs. """

  def __init__(self):
    super(JobDict, self).__init__()
    # List of subjobs (as in subdirectories). 
    super(JobDict, self).__setattr__("children", {})
    # This particular job. 
    super(JobDict, self).__setattr__("params", {})
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
    if value is not None and not hasattr(value, "__call__"):
      raise ValueError("job.functional should be either None(no job) or a callable.")
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
     if self.parent is None: return "/"
     string = None
     for key, item in self.parent.children.iteritems():
       if id(item) == id(self):
         string = self.parent.name + key
         break
     if string is None: raise RuntimeError("Could not determine the name of the dictionary.")
     return string + '/'

  def __getitem__(self, index): 
    """ Returns job description from the dictionary.

        If the job does not exist, will create it.
    """
    from re import split
    from os.path import normpath

    index = normpath(index)
    if index == "" or index is None or index == ".": return self
    if index[0] == "/": return self.root[index[1:]]

    result = self
    names = split(r"(?<!\\)/", index)
    for i, name in enumerate(names):
      if name == "..":
        if result.parent is None: raise KeyError("Cannot go below root level.")
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
      if id(self) == id(deletee): raise KeyError("Will not commit suicide.")
      parent = self.parent
      while parent is not None: 
        if id(parent) == id(deletee): raise KeyError("Will not go Oedipus on you.")
        parent = parent.parent

    parent = self[index+"/.."]
    name = relpath(index, index+"/..")
    if name in parent.children:
      if id(self) == id(parent.children[name]): raise KeyError("Will not delete self.")
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
      if parentpath not in self:
        raise KeyError('Could not find parent job {0}.'.format(parentpath))
      mother = self[parentpath]
      parent = self.parent
      while parent is not None:
        if parent is mother: raise KeyError('Will not set parent job of current job.')
    if len(childpath) == 0 or childpath == '.': raise KeyError('Will not set current directory.')
    if childpath == '..': raise KeyError('Will not set parent directory.')

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
      while result.parent is not None: result = result.parent
      return result / index[1:]

    names = split(r"(?<!\\)/", index) 
    result = self
    for name in names:
      if name == "..":
        if result.parent is None:
          raise RuntimeError('Cannot descend below root.')
        result = result.parent
        continue
      elif name not in result.children:
        result.children[name] = JobDict()
        result.children[name].parent = result
      result = result.children[name]
    return result

  @property
  def is_job(self):
    """ True if functional is not None. """
    return self.functional is not None

  def subjobs(self):
    """ Iterator over children jobs. """
    return sorted(self.children.iterkeys())
    
  def compute(self, **kwargs):
    """ Performs calculations over job list. """  

    if not self.is_job: return None
    kwargs.update(self.params)
    return self.functional.__call__(**kwargs)

  def update(self, other, merge=False):
    """ Updates job and tree with other.
    
        :Parameters: 
           other : `JobDict`
             An other job tree from which to update.
           merge : bool
             If false (default), then actual jobs in ``other`` completely
             overwrite actual jobs in ``self``. If False, then ``params`` in
             ``self`` is updated with ``params`` in ``other`` if either one is
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
      self.params = other.params
      self.functional = other.functional
    else:
      if not (self.is_job or other.is_job): return
      self.params.update(other.params)
      if other.functional is not None: self.functional = other.functional

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
    if name in self.params: return self.params.pop(name)
    raise AttributeError("Unknown job attribute " + name + ".")

  def __getattr__(self, name):
    """ Returns job attribute. """
    if name in self.params: return self.params[name]
    raise AttributeError("Unknown job attribute " + name + ".")

  def __setattr__(self, name, value):
    """ Sets job attribute. """
    from pickle import dumps
    if name in self.params:
      try: dumps(value)
      except Exception as e:
        raise ValueError("Could not pickle job-parameter. Caught error:\n{0}".format(e))
      else: self.params[name] = value
    else: super(JobDict, self).__setattr__(name, value)

  def __dir__(self):
    from itertools import chain
    result = chain([u for u in self.__dict__ if u[0] != '_'], \
                   [u for u in dir(self.__class__) if u[0] != '_'], \
                   [u for u in self.params.iterkeys() if u[0] != '_'])
    return list(set(result))

  def __getstate__(self):
    d = self.__dict__.copy()
    params = d.pop("params")
    return d, params, "This is a JobDict pickle. Grep me!"
  def __setstate__(self, args):
    super(JobDict, self).__setattr__("params", args[1])
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
  __iter__ = iterkeys
  """ Iterator over keys. """


  @property
  def nbjobs(self):
    """ Returns the number of jobs in tree. """
    return len([0 for j, o in self.iteritems() if not j.is_tagged])

  @property 
  def root(self): 
    """ Returns root dictionary. """
    result = self
    while result.parent is not None: result = result.parent
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
        params. However, functional and params values should the same
        object as self. The sub-branches of the returned dictionary are shallow
        copies of the sub-branches of self. In other words, the functional and
        refences in params dictionary are in common between result and self,
        but nothing else.

        The returned dictionary does not have a parent!
    """
    from copy import copy
    # new job-dictionary.
    result = JobDict()
    result._functional = self._functional
    result.params   = self.params.copy()
    result.parent     = None
    for name, value in self.children.items():
      result.children[name] = copy(value)
      result.children[name].parent = result
    attrs = self.__dict__.copy()
    attrs.pop('params')
    attrs.pop('parent')
    attrs.pop('children')
    attrs.pop('_functional')
    result.__dict__.update(attrs)
    return result


class SuperCall(object):
  """ Obviates issues when using a "super" functional.

      Since functionals of a jobdictionary are deepcopied, the following line
      will not result in calling the next class in the __mro__.
  
      >>> jobdict.functional = super(Functional, functional)

      Indeed, this line will first call the __getitem__, __setitem__ (or
      __deepcopy__) of the super object. In general, this means we end-up with
      ``jobdict.function == functional``.

      This class obviates this difficulty.

      >>> jobdict.functional = SuperCall(Functional, functional)
  """
  def __init__(self, class_, object_):
    object.__init__(self)
    self.__dict__['_class'] = class_
    self.__dict__['_object'] = object_
  def __call__(self, *args, **kwargs):
    return super(self._class, self._object).__call__(*args, **kwargs)
  def __getattr__(self, name):
    try: return getattr(super(self._class, self._object), name)
    except: return getattr(self._object, name)
  def __setattr__(self, name, value): setattr(self._object, name, value)
  def __getstate__(self): return self._class, self._object
  def __setstate__(self, args):
    self.__dict__['_class'] = args[0]
    self.__dict__['_object'] = args[1]
  def __dir__(self): return dir(super(self._class, self._object))
  
    

