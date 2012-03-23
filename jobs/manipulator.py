""" Classes to manipulate job-dictionaries. """
__docformat__ = "restructuredtext en"
__all__ = ['JobParams']

from .extract import AbstractMassExtract

class JobParams(AbstractMassExtract):
  """ Get and sets job parameters for a job-dictionary. """
  def __init__(self, jobdict=None, only_existing=True, **kwargs):
    """ Initializes job-parameters.

        :param jobdict :
            :py:class:`JobDict <lada.jobs.jobdict.JobDict>` instance for which
            to get/set parameters. If None, will look for ipython's
            current_jobdict.
        :param bool only_existing:
            If true (default), then only existing parameters can be modified:
            non-existing parameters will not be added.
        :param kwargs:
            Variable length keyword argument passed on to
            :py:meth:`AbstractMassExtract.__init__`.
    """
    from .. import jobparams_only_existing

    self.__dict__["_jobdict"] = jobdict
    self.__dict__['only_existing'] = jobparams_only_existing\
                                       if only_existing is None\
                                       else only_existing
    """ Only modifies parameter which already exist. """

    if 'dynamic' not in kwargs: kwargs['dynamic'] = True
    super(JobParams, self).__init__(**kwargs)

  @property
  def jobdict(self):
    from lada.jobs import load
    from lada import is_interactive
    if self._jobdict is None:
      if self._rootpath is None: 
        if is_interactive:
          from lada import interactive
          if interactive.jobdict is None:
            print "No current job-dictionary."
            return
          return interactive.jobdict
        else: raise RuntimeError('No jobdictionary.')
      else: self._jobdict = load(self.rootpath, timeout=30)
    return self._jobdict.root


  @property
  def addattr(self):
    """ Returns manipulator with ability to *add new* attributes. """
    return self.shallow_copy(only_existing=False)
    
  @property
  def onoff(self):
    """ Dictionary with calculations which will run.

	Whereas other properties only report untagged jobs, this will report
        both. Effectively checks wether a job is tagged or not. Calculations which 
    """
    result = {}
    for name, job in self.iteritems():
      result[name] = "off" if job.is_tagged else "on"
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return result

  @onoff.setter
  def onoff(self, value):
    """ Dictionary with tagged and untagged jobs.

        Whereas other properties only report untagged jobs, this will report
        both.
    """
    if hasattr(value, 'iteritems'):
      for key, value in value.iteritems():
        try: job = self[key]
        except: continue
        else: job.onoff = value
    elif value == "on" or value == True:
      for name, job in self.iteritems(): job.untag()
    elif value == "off" or value == False:
      for name, job in self.iteritems(): job.tag()

  @property
  def extractors(self):
    """ Returns dictionary of extrators. """
    from .forwarding_dict import ForwardingDict
    result = self.dicttype()
    for k, j in self.iteritems(): result[k] = j
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return ForwardingDict( dictionary=result, naked_end=self.naked_end, \
                           only_existing=self.only_existing, readonly=False)
    

  @property
  def view(self):
    """ A regex pattern which the name of extracted jobs should match.

        If None, then no match required. Should be a string, not an re object.
    """
    return self.jobdict.name if self._view is None else self._view
  @view.setter
  def view(self, value): self._view = value

  def __iter_alljobs__(self):
    """ Loops through all correct jobs. """
    for name, job in self.jobdict.iteritems(): yield job.name, job

  def __getattr__(self, name): 
    """ Returns extracted values. """
    from .forwarding_dict import ForwardingDict
    result = self.dicttype()
    for key, value in self.iteritems():
      if value.is_tagged: continue
      try: result[key] = getattr(value, name)
      except: result.pop(key, None)
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    if len(result) == 0: 
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )
    return ForwardingDict( dictionary=result, naked_end=self.naked_end, \
                           only_existing=self.only_existing, readonly=False)

  def __setattr__(self, name, value):
    """ Returns dictionary with job parameters for each job. """
    from re import match
    # initialization not done yet.
    if "only_existing" not in self.__dict__: super(JobParams, self).__setattr__(name, value)
    # some cached attributes.
    if match("_cached_attr\S+", name): super(JobParams, self).__setattr__(name, value)
    # Look for other attriubtes in current instance.
    try: super(JobParams, self).__getattribute__(name)
    except AttributeError: pass
    else:
      super(JobParams, self).__setattr__(name, value)
      return 

    found = False
    for jobname, job in self.iteritems():
      if job.is_tagged: continue
      if hasattr(job, name):
        setattr(job, name, value)
        found = True
      elif not self.only_existing: 
        job.jobparams[name] = value
        found = True
    if not found:
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )

  def __delattr__(self, name):
    try: super(JobParams, self).__getattribute__(name)
    except AttributeError: pass
    else:
      super(JobParams, self).__delattr__(name, value)
      return

    found = False
    for jobname, job in self.iteritems():
      if job.is_tagged: continue
      if hasattr(job, name):
        delattr(job, name)
        found = True
    if not found:
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )

  @property
  def _attributes(self):
    """ Attributes which already exist. """
    result = set()
    for name, job in self.iteritems():
      if not job.is_tagged: result |= set([u for u in dir(job) if u[0] != '_'])
    return result

  def __setitem__(self, name, jobdict):
    """ Modifies/creates item in jobdictionary.
    
        :param str name: 
          Name of item to modify.
        :param jobdict:
          :py:class:`JobDict lada.jobs.jobdict.JobDict` or
          :py:class:`JobParams` with which to set/modify item.
          In the latter case, it should point to a single entry in the
          jobdictionary. Eg no wildcards.

        This function provides the ability to extend a jobdictionary with other jobs.

        >>> jobparams['newjob'] = jobparams['oldjob'] 
        >>> jobparams['newjob'] = some_job_dictionary

        In  both cases above, the left-hand-side cannot be a wildcard.
        Similarly, in the first case above, the right hand side should also point to
        a valid job, not a sequence of jobs:

        >>> jobparams['newjobs'] = jobparams['*/oldjobs']
        raises KeyError
        >>> jobparams['*/newjobs'] = jobparams['oldjobs']
        No Error! creates a job called '*/newjobs'

        .. warning:: The right-hand-side is *always* deep-copied_.
          .. _deep-copied:: http://docs.python.org/library/copy.html
    """
    from lada import is_interactive
    from copy import deepcopy

    if isinstance(jobdict, JobParams): 
      try: jobdict = jobdict.jobdict[jobdict.view]
      except KeyError: raise KeyError("{0} is not an actual job".format(jobdict.view))
    if name in self.jobdict and is_interactive:
      print "Modifying existing job {0}.".format(name)
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("Modifying existing job {0}.\nIs this OK? [y/n] ".format(name))
      if a == 'n':
        print "Aborting."
        return
    self.jobdict[name] = deepcopy(jobdict)
 
  def __delitem__(self, name):
    """ Deletes items from jobdictionary. """
    from lada import is_interactive
    if is_interactive:
      print "Deleting the following jobs:"
      for key in self[name].keys(): print key
      a = ''
      while a != 'n' and a != 'y':
        a = raw_input('Ok? [y/n] ')
      if a == 'n':
        print "Aborting."
        return
    for key in self[name].keys(): del self.jobdict.root[key]

  def concatenate(self, jobdict):
    """ Updates content of current jobdictionary with that of another.

        :param jobdict:
          :py:class:`JobDict` instance, or :py:class:`JobParams` instance with
          which to update the current job-dictionary.
        
        Update means that jobs and jobparameters will be overwritten with those
        from the input. Jobs in the input which are not in the current
        job-dictionary will be overwritten. If `jobdict` is a
        :py:class:`JobDict` instance, it is possible to use wildcards in order
        to select those jobs of interests.

        .. warning: New jobs are always added at the root of the job-dictionary.
          Make sure the jobs bear the names you want.
    """
    from .jobdict import JobDict
    from lada import is_interactive
    keys = jobdict.keys()
    if is_interactive:
      if len(keys) == 0: 
        print "Empty input jobdictionary. Aborting."
        return
      add = [k for k in keys if k in self]
      if len(add) > 0:
        print "Adding the following jobdictionaries:"
        for key in add: print key
      update = [k for k in keys if k in self]
      if len(update) > 0:
        print "Updating the following jobdictionaries:"
        for key in update: print key
      a = ''
      while a != 'n' and a != 'y':
        a = raw_input("Is the above OK? [n/y] ")
      if a == 'n':
        print "Aborting."
        return
    rootadd = jobdict
    if isinstance(rootadd, JobParams):
      rootadd = JobDict()
      for key in keys:
        job = rootadd / key 
        rootadd[key] = jobdict.jobdict[key]
    
    self.jobdict.root.update(rootadd)


