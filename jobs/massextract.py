""" Classes to manipulate output from jobdictionaries. """
__docformat__ = "restructuredtext en"
__all__ = ['MassExtract']
from .extract import AbstractMassExtract

class MassExtract(AbstractMassExtract): 
  """ Propagates extraction methods from different jobs. 
  
      Collects extractors across all jobs (for which job.functional.Extract
      exist). The results are presented as attributes of an instance of
      MassExtract, and arranged as directory where the key is the name of the
      job and the value obtained from an instance of that job's Extract. This
      class is set-up to fail silently, and hence is of limited use for
      diagnosis.
  """

  def __init__(self, path=None, **kwargs):
    """ Initializes extraction object. 
 
        :param str path:
            Pickled jobdictionary for which to extract stuff. If None, will
            attempt to use the current jobdictionary.
        :param kwargs:
            Variable length keyword argument passed on to
            :py:meth:`AbstractMassExtract.__init__`.
    """
    from ..misc import RelativePath
    self.__dict__["_jobdict"] = None
    super(MassExtract, self).__init__(path=path, **kwargs)

  @property
  def view(self):
    """ A regex pattern which the name of extracted jobs should match.

        If None, then no match required. Should be a string, not an re object.
    """
    return self._view if self._view is not None else self.jobdict.name
  @view.setter
  def view(self, value): self._view = value

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

  def __iter_alljobs__(self):
    """ Generator to go through all relevant jobs.  
    
        :return: (name, extractor), where name is the name of the job, and
          extractor an extraction object.
    """
    from os.path import join, dirname
    
    for name, job in self.jobdict.iteritems():
      if job.is_tagged: continue
      try: extract = job.functional.Extract(join(dirname(self.rootpath), name))
      except: pass 
      else: yield job.name, extract
