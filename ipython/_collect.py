""" Subpackage for collect mechanism """
from ..jobs import AbstractMassExtract
__all__ = ['Collect']
class Collect(AbstractMassExtract):
  """ Mass extraction with varying position argument. 
  
      By adjusting ``self.position``, which jobs to collect can be adjusted.
  """
  def __init__(self, comm=None, _view=None, naked_end=True):
    """ Initializes a Collect instance. """
    from IPython.ipapi import get as get_ip_handle
    super(Collect, self).__init__(None, comm=comm)

    self.ip = get_ip_handle()
    """ Gets current handle. """
    self._view = _view
    """ Whether this is a view, or the original collect. 
    
        Should be None if not a view. In that case, the Collect instance will
        look for the position of the current_jobdict. Otherwise, this should be
        a string which forms the beginning of the job-dictionaries to collect.
    """
    self.naked_end = naked_end
    """ If True, a value, rathe than a dict, is returned if at end of branch. """

  def walk_through(self):
    """ Generator to go through all relevant jobs.  
    
        :return: (name, extractor), where name is the name of the job, and
          extractor an extraction object.
    """
    from os.path import exists, join
    
    if not self.root: return
    jobdict = self.jobdict
    for job, name in jobdict.walk_through():
      if job.is_tagged: continue
      if not hasattr(job.functional, "Extract"): continue
      try: extract = job.functional.Extract(join(self.root, name), comm = self.comm)
      except: pass
      else: yield name, extract

  @property 
  def jobdict(self):
    """ Returns root of current dictionary. """
    if "current_jobdict" not in self.ip.user_ns: 
      print "No current job-dictionary to collect from."
      return
    return self.ip.user_ns["current_jobdict"].root

  @property 
  def position(self):
    """ Returns current position in dictionary. """
    if "current_jobdict" not in self.ip.user_ns: 
      print "No current job-dictionary to collect from."
      return
    return (self._view if self._view != None else self.ip.user_ns["current_jobdict"].name)[1:]

  @property 
  def root(self):
    """ Returns directory name of current dictionary. """
    from os.path import dirname
    if "current_jobdict_path" not in self.ip.user_ns: 
      print "No known path for current job-dictionary."
      print "Don't know where to look for results."
      print "Please use %saveto magic function."
      return
    return dirname(self.ip.user_ns["current_jobdict_path"])

  def __getattr__(self, name): 
    """ Returns extracted values. """
    from os.path import dirname
    from re import compile
    if name == "_cached_extractors" or name == "_cached_properties": 
      raise AttributeError("Unknown attribute {0}.".format(name))
    if name in self._properties(): 
      result = {}
      position = compile("^" + self.position)
      for key, value in self._extractors().items():
        if position.match(key) == None: continue
        try: result[key] = getattr(value, name)
        except: result.pop(key, None)
      if self.naked_end and len(result.keys()) == 1: return result[result.keys()[0]] 
      return result
    raise AttributeError("Unknown attribute {0}.".format(name))

  def __getitem__(self, name):
    """ Returns a view of the current job-dictionary. """
    from os.path import normpath, join
    if name[0] == '/': return self.__class__(comm=self.comm, _view=name)
    path = normpath(join('/', join(self.position, name)))
    return self.__class__(comm=self.comm, _view=path)

  @property
  def children(self):
    """ Iterates through all sub-jobs. """
    from os.path import join, normpath
    try: jobdict = self.jobdict["/" + self.position]
    except: return
    for name in jobdict.children.keys():
      path = normpath(join(join("/", self.position), name))
      yield self.__class__(comm=self.comm, _view=path)

  def grep(self, regex, flags=0):
    """ Yields views for children with fullnames matching the regex.
    
        :Param regex: The matching regular expression.
        :type regex: str

        The match is successful if the regex is matched using python's
        `re.search <http://docs.python.org/library/re.html#re.search>`_ method.

        Only the outermost view of each math is given. In other words, if a
        view is yielded, its subviews will not be yielded.
    """
    from re import compile
    reg = compile(regex, flags)
    if reg.search(self.position) != None:
      yield self
      return
    for child in self.children:
      if reg.search(child.position) != None: yield child
      else: # goes to next level.
        for grandchild in child.grep(regex, flags): yield grandchild
