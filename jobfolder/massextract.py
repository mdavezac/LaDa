""" Classes to manipulate output from jobfolderionaries. """
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
            Pickled job-folder for which to extract stuff. If None, will
            attempt to use the current job-folder.
        :param kwargs:
            Variable length keyword argument passed on to
            :py:meth:`AbstractMassExtract.__init__`.
    """
    self.__dict__["_jobfolder"] = None
    super(MassExtract, self).__init__(path=path, **kwargs)

  @property
  def view(self):
    """ A regex pattern which the name of extracted jobs should match.

        If None, then no match required. Should be a string, not an re object.
    """
    return self._view if self._view is not None else self.jobfolder.name
  @view.setter
  def view(self, value): self._view = value

  @property
  def jobfolder(self):
    from . import load
    from .. import is_interactive
    if self._jobfolder is None:
      if self._rootpath is None: 
        if is_interactive:
          from .. import interactive
          if interactive.jobfolder is None:
            print "No current job-dictionary."
            return
          return interactive.jobfolder
        else: raise RuntimeError('No job-folder.')
      else: self._jobfolder = load(self.rootpath, timeout=30)
    return self._jobfolder.root

  @property
  def rootpath(self):
    from .. import is_interactive
    if self._jobfolder is None and self._rootpath is None and is_interactive:
      from .. import interactive
      if interactive.jobfolder_path is None:
        print "No current path to job-dictionary."
        return
      return interactive.jobfolder_path
    return super(MassExtract, self).rootpath
    

  def __iter_alljobs__(self):
    """ Generator to go through all relevant jobs.  
    
        :return: (name, extractor), where name is the name of the job, and
          extractor an extraction object.
    """
    from os.path import join, dirname
    
    for name, job in self.jobfolder.iteritems():
      if job.is_tagged: continue
      try: extract = job.functional.Extract(join(dirname(self.rootpath), name))
      except: pass 
      else: yield job.name, extract
