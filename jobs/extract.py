""" Classes to manipulate output from jobdictionaries. """
__docformat__ = "restructuredtext en"
__all__ = ['AbstractMassExtract', 'MassExtract', 'AbstractMassExtractDirectories']

from abc import ABCMeta, abstractmethod
from .forwarding_dict import ForwardingDict

class AbstractMassExtract(object): 
  """ Propagates extraction methods from different jobs. """
  __metaclass__ = ABCMeta

  def __init__(self, view=None, excludes=None, dynamic=False, ordered=True, **kwargs):
    """ Initializes extraction object. 

        :Parameters:
          view : str or None
            Pattern which the job names must match to be included in the
            extraction.
          excludes : list of str or None
            List of patterns which the job names must *not* match to be
            included in the extraction.
          dynamic : boolean
            If true, chooses a slower but more dynamic caching method. Only
            necessary for ipython shell. 
          ordered : boolean
            If true, uses OrderedDict rather than conventional dict.

        :Kwarg naked_end: True if should return value rather than dict when only one item.
        :Kwarg unix_re: converts regex patterns from unix-like expression.
    """
    from . import default_params as default
    from ..opt import OrderedDict

    object.__init__(self)

    self.naked_end = kwargs.pop('naked_end', default.naked_end)
    """ If True and dict to return contains only one item, returns value itself. """
    self.view = view
    """ The pattern which job-names should match. """
    self.unix_re = kwargs.pop('unix_re', default.unix_re)
    """ If True, then all regex matching is done using unix-command-line patterns. """
    self.excludes = excludes
    assert len(kwargs) == 0, ValueError("Unkwnown keyword arguments:{0}.".format(kwargs.keys()))
    self._cached_extractors = None
    """ List of extration objects. """
    self.dynamic = dynamic
    """ If True chooses a slower but more dynamic caching method. """
    self.dicttype = OrderedDict if ordered else dict

  def uncache(self): 
    """ Uncache values. """
    self._cached_extractors = None

  @property 
  def excludes(self):
    """ Pattern or List of patterns to ignore. or None.

        ``self.unix_re`` determines whether these are unix-command-line like
        patterns or true python regex.
    """ 
    try: return self._excludes 
    except AttributeError: return None
  @excludes.setter
  def excludes(self, value):
    if isinstance(value, str): self._excludes = [value]
    else: self._excludes = value

  def avoid(self, excludes):
    """ Returns a new MassExtract object with further exclusions. 

        :Param excludes: Pattern or patterns to exclude from output.
        :type excludes: str or list of str or None 
          
        The goal of this function is to work as an *anti* operator [], i.e. by
        excluding from the output anything that matches the patterns, rather
        including only those which match the pattern.
        This is strickly equivalent to:

        >>> other = massextract.copy(excludes=excludes)
        >>> other.excludes.extend(massextract.excludes)

        and then doing calculations with ``other``. The advantage is that it
        can all be done on one line.

        If the ``excludes`` argument is None or an empty list, then the
        returned object will not exlude anything.
    """ 
    if excludes == None or len(excludes) == 0: return self.copy(excludes=None)
    result = self.copy(excludes=excludes)
    if self.excludes != None: result.excludes.extend(self.excludes)
    return result

  def __iter__(self):
    """ Iterates through all job names. """
    for name, job in self._regex_extractors(): yield name

  def iteritems(self):
    """ Iterates through all extraction objects and names. """
    for name, job in self._regex_extractors(): yield name, job
  def items(self):
    """ Iterates through all extraction objects and names. """
    return [(name, job) for name, job in self.iteritems()]
    
  def itervalues(self):
    """ Iterates through all extraction objects. """
    for name, job in self._regex_extractors(): yield job
  def values(self):
    """ Iterates through all extraction objects. """
    return [job for job in self.itervalues()]

  def iterkeys(self):
    """ Iterates through all extraction objects. """
    for name, job in self._regex_extractors(): yield name
  def keys(self):
    """ Iterates through all extraction objects. """
    return [name for name in self.iterkeys()]
  
  def __iter__(self):
    """ Iterates through all job names. """
    for name, job in self.iteritems(): yield name
  def __len__(self): 
    """ Returns length of output dictionary. """
    return len(self.keys())

  def _regex_pattern(self, pattern, flags=0):
    """ Returns a regular expression. """
    from re import compile
    from ..opt import convert_from_unix_re
    return compile(pattern, flags) if not self.unix_re\
           else convert_from_unix_re(pattern)

  @abstractmethod
  def __iter_alljobs__(self):
    """ Generator to go through all relevant jobs. 
    
        :return: (name, extractor), where name is the name of the job, and
          extractor an extraction object.
    """
    pass

  @property
  def view(self):
    """ A regex pattern which the name of extracted jobs should match.

        If None, then no match required. Should be a string, not an re object.
    """
    if self._view == None: return ""
    return self._view
  @view.setter
  def view(self, value): self._view = value

  @property
  def _extractors(self):
    """ Goes through all jobs and collects Extract if available. """
    if self.dynamic:
      if self._cached_extractors == None: self._cached_extractors = self.dicttype()
      result = self.dicttype()
      for name, extract in self.__iter_alljobs__():
        if name not in self._cached_extractors: self._cached_extractors[name] = extract
        result[name] = self._cached_extractors[name]
      return result
    else:
      if self._cached_extractors != None: return self._cached_extractors
      result = self.dicttype()
      for name, extract in self.__iter_alljobs__(): result[name] = extract
      self._cached_extractors = result
      return result

  def _regex_extractors(self):
    """ Loops through jobs in this view. """
    if self._view == "": 
      for key, value in self._extractors.iteritems(): yield key, value
      return

    regex = self._regex_pattern(self.view)
    if self.excludes != None: excludes = [self._regex_pattern(u) for u in self.excludes]
    for key, value in self._extractors.iteritems():
      if regex.match(key) == None: continue
      if self.excludes != None and any(u.match(key) != None for u in excludes): continue
      yield key, value

  @property
  def _attributes(self): 
    """ Returns __dir__ special to the extraction itself. """
    results = set([])
    for key, value in self.iteritems():
      results |= set([u for u in dir(value) if u[0] != '_'])
    return results

  def __dir__(self): 
    from itertools import chain
    results = chain( [u for u in self.__dict__ if u[0] != '_'], \
                     [u for u in dir(self.__class__) if u[0] != '_'], \
                     self._attributes )
    return list(set(results))

  def __getattr__(self, name): 
    """ Returns extracted values. """
    assert name in self._attributes, AttributeError("Unknown attribute {0}.".format(name))

    result = self.dicttype()
    for key, value in self.iteritems():
      try: result[key] = getattr(value, name)
      except: result.pop(key, None)
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return ForwardingDict(result, naked_end=self.naked_end)

  def __getitem__(self, name):
    """ Returns a view of the current job-dictionary. """
    from os.path import normpath, join
    if name[0] == '/': return self.copy(view=name)
    path = normpath(join('/', join(self.view, name)))
    return self.copy(view=path)

  @property
  def children(self):
    """ next set of minimal regex. """
    from os.path import join, normpath
    regex = self._regex_pattern(self.view)

    jobs = self.keys()
    if len(self.keys()) < 2: return 
    children = set()
    if len(self.view) == 0 or self.view == '/':
      for name in self.iterkeys():
        children.add(name[:1+name[1:].find('/')])
    else:
      for name in self.iterkeys():
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
          yield_match : bool
             If True, will yield a two tuple, where the second item is the
             match object.
             If False, only the view is yielded.
             This option is not available (or meaningfull) if ``self.unix_re``
             is True.

        The match is successful if the regex is matched using python's
        `re.search`__ method.

        .. __:  http://docs.python.org/library/re.html#re.search

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
        only, eg without blocking mpi calls. `solo` returns an
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
    result = self.__class__()
    result.__dict__.update(self.__dict__)
    return result

  def copy(self, **kwargs):
    """ Returns a shallow copy. 
    
        :Param kwargs:  Any keyword attribute will modify the corresponding
          attribute of the copy.
    """
    from copy import copy
    result = copy(self)
    for key, value in kwargs.iteritems(): setattr(result, key, value)
    return result


  @property
  def jobs(self):
    """ Deprecated. Use keys and iterkeys instead. """
    from warnings import warn
    warn('jobs property is deprecated. Please use keys and or iterkeys instead.', DeprecationWarning)
    return self.keys()



class MassExtract(AbstractMassExtract): 
  """ Propagates extraction methods from different jobs. 
  
      Collects extractors across all jobs (for which job.functional.Extract
      exist). The results are presented as attributes of an instance of
      MassExtract, and arranged as directory where the key is the name of the
      job and the value obtained from an instance of that job's Extract. This
      class is set-up to fail silently, and hence is of limited use for
      diagnosis.
  """

  def __init__(self, path=None, comm=None, **kwargs):
    """ Initializes extraction object. 
 
        :Parameters:
          path : str or None
            Pickled jobdictioanary for which to extract stuff. If None, will
            attempt to use the current jobdictionary.
          comm : boost.mpi.communicator or None
            Optional communicator. How communicators are used will depend on
            each calculation's extractor.
          kwargs : dict
            Variable length keyword argument passed on to `AbstractMassExtract`.

        :kwarg view: Pattern to match to job names.
        :kwarg excludes: List of patterns which job-names should not match.
        :kwarg naked_end: True if should return value rather than dict when only one item.
        :kwarg unix_re: converts regex patterns from unix-like expression.
    """
    AbstractMassExtract.__init__(self, **kwargs)

    from os.path import isdir, isfile, exists, dirname, abspath

    self.rootdir = path # Path to the job dictionary.
    self.comm = comm

  @property
  def comm(self):
    """ MPI Communicator, or None for serial. 

        This property is intended to synchronize communicator over all
        extractor objects. How MPI is done will depend on individual
        extractors. Note that extractors are initialized with communicators
        only if they accept a ``comm`` keyword. Communicators are set only if
        an extractor contains a ``comm`` attribute.
    """
    return self._comm
  @comm.setter
  def comm(self, value):
    if self._cached_extractors != None:
      for e in self._cached_extractors.values(): e.comm = value
    self._comm = value
  @comm.deleter
  def comm(self):
    if self._cached_extractors != None:
      for e in self._cached_extractors.values():
        if hasattr(e, "comm"): e.comm = None
    self._comm = None

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

  def __iter_alljobs__(self):
    """ Generator to go through all relevant jobs.  
    
        :return: (name, extractor), where name is the name of the job, and
          extractor an extraction object.
    """
    from os.path import join
    
    for name, job in self.jobdict.iteritems():
      if job.is_tagged: continue
      try: extract = job.functional.Extract(join(self.rootdir, name), comm = self.comm)
      except: pass 
      else: yield job.name, extract

class AbstractMassExtractDirectories(AbstractMassExtract):
  """ Propagates extractors from all subdirectories.
  
      Trolls through all subdirectories for calculations with given extraction
      files, and organises results as a dictionary where keys are the name of
      the diretory.

      An class derived from this one should make sure that:
      
      - `Extract` is not none.
      - `__is_calc_dir__ ` is correctly defined. 
  """
  def __init__(self, path = '.', Extract = None, comm = None, **kwargs):
    """ Initializes AbstractMassExtractDirectories.
    
    
        :Parameters:
          path : str 
            Root directory for which to investigate all subdirectories.
            If None, uses current working directory.
          Extract
            Extraction class to use within each calculation. 
          comm : boost.mpi.communicator or None
            Optional communicator. How communicators are used will depend on
            each calculation's extractor.
          kwargs : dict
            Keyword parameters passed on to AbstractMassExtract.

        :kwarg naked_end: True if should return value rather than dict when only one item.
        :kwarg unix_re: converts regex patterns from unix-like expression.
    """
    from os.path import exists, isdir
    from ..opt import RelativeDirectory

    # this will throw on unknown kwargs arguments.
    AbstractMassExtract.__init__(self,**kwargs)

    self.Extract = Extract
    """ Extraction class to use. """

    self._rootdir = RelativeDirectory(path, hook=self.uncache)
    """ Root of the directory-tree to trawl for OUTCARs. """
    
    # mpi communicator is a property.
    self.comm = comm

  @property
  def comm(self):
    """ MPI Communicator, or None for serial. 

        This property is intended to synchronize communicator over all
        extractor objects. How MPI is done will depend on individual
        extractors. Note that extractors are initialized with communicators
        only if they accept a ``comm`` keyword. Communicators are set only if
        an extractor contains a ``comm`` attribute.
    """
    return self._comm
  @comm.setter
  def comm(self, value):
    if self._cached_extractors != None:
      for e in self._cached_extractors.values(): e.comm = value
    self._comm = value
  @comm.deleter
  def comm(self):
    if self._cached_extractors != None:
      for e in self._cached_extractors.values():
        if hasattr(e, "comm"): e.comm = None
    self._comm = None

  @property
  def rootdir(self): 
    """ Root of the directory-tree to trawl for OUTCARs. """
    return self._rootdir.path
  @rootdir.setter
  def rootdir(self, value): self._rootdir.path = value

  def __iter_alljobs__(self):
    """ Goes through all directories with a contcar. """
    from os import walk, getcwd
    from os.path import abspath, relpath, abspath, join

    for dirpath, dirnames, filenames in walk(self.rootdir, topdown=True, followlinks=True):
      if not self.__is_calc_dir__(dirpath, dirnames, filenames): continue

      try: result = self.Extract(join(self.rootdir, dirpath), comm = self.comm)
      except TypeError: # no comm keyword.  
        try: result = self.Extract(join(self.rootdir, dirpath))
        except: continue
      except: continue

      result.OUTCAR = self.OUTCAR
      yield join('/', relpath(dirpath, self.rootdir)), result

  @property
  def _attributes(self): 
    """ Returns __dir__ set special to the extraction itself. """
    return set([u for u in dir(self.Extract()) if u[0] != '_'])
  
  @abstractmethod
  def __is_calc_dir__(self, dirpath, dirnames, filenames):
    """ Returns true this directory contains a calculation. """
    pass

  def __copy__(self):
    """ Returns a shallow copy of this object. """
    from copy import copy
    result = self.__class__(self.rootdir)
    for k, v in self.__dict__.items():
      if k != '_rootdir': result.__dict__[k] = copy(v)
    return result
