""" Subpackage for collect mechanism """
from ..jobs import AbstractMassExtract
__all__ = ['Collect']
class Collect(AbstractMassExtract):
  """ Mass extraction with varying position argument. 
  
      By adjusting ``self.position``, which jobs to collect can be adjusted.
  """
  def __init__(self, comm=None, _view=None):
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

  def walk_through(self):
    """ Generator to go through all relevant jobs.  
    
        :return: (name, extractor), where name is the name of the job, and
          extractor an extraction object.
    """
    from os.path import exists, join
    
    for job, name in self.jobdict.walk_through():
      if job.is_tagged: continue
      if not hasattr(job.functional, "Extract"): continue
      if not exists(join(self.root, name)): continue
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
    if self._view != None: return self._view[1:]
    if "current_jobdict" not in self.ip.user_ns: 
      print "No current job-dictionary to collect from."
      return
    return self.ip.user_ns["current_jobdict"].name[1:]

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
    if name == "_cached_extractors" or name == "_cached_properties": 
      raise AttributeError("Unknown attribute {0}.".format(name))
    if name in self._properties(): 
      result = {}
      position = self.position
      for key, value in self._extractors().items():
        if len(position) > 0:
          if position != key[:len(position)]: continue
          if len(position) < len(key) and position[-1] != '/': continue
        try: result[key] = getattr(value, name)
        except: result.pop(key, None)
      return result if len(result.keys()) > 1 else result[result.keys()[0]]
    raise AttributeError("Unknown attribute {0}.".format(name))

  def __getitem__(self, name):
    """ Returns a view of the current job-dictionary. """
    if name[0] == '/': return self.__class__(comm=self.comm, _view=name)
    return self.__class__(comm=self.comm, _view="/"+self.position+name)

  @property
  def children(self):
    """ Iterates through all sub-jobs. """
    if "/" + self.position not in self.jobdict.children: return
    for name in self.jobdict["/"+self.position].children.keys():
      yield self.__class__(comm=self.comm, _view="/"+self.position + name)
