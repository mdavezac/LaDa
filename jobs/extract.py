""" Classes to manipulate output from jobfolderionaries. """
__docformat__ = "restructuredtext en"
__all__ = ['AbstractMassExtract', 'MassExtract', 'AbstractMassExtractDirectories']

from abc import ABCMeta, abstractmethod

class AbstractMassExtract(object): 
  """ Propagates extraction methods from different jobs. """
  __metaclass__ = ABCMeta

  def __init__(self, path=None, view=None, excludes=None, dynamic=False, ordered=True, 
               naked_end=None, unix_re=True):
    """ Initializes extraction object. 

        :param str path:
            Root directory for which to investigate all subdirectories.
            If None, uses current working directory.
        :param str view:
            Pattern which the job names must match to be included in the
            extraction. Ignored if None.
        :para excludes:
            List of patterns which the job names must *not* match to be
            included in the extraction. Ignored if None.
        :param bool dynamic:
            If true, chooses a slower but more dynamic caching method. Only
            usefull for ipython shell. 
        :param ordered : boolean
            If true, uses OrderedDict rather than conventional dict.
        :param bool naked_end:
            True if should return value rather than dict when only one item.
        :param bool unix_re: 
            Converts regex patterns from unix-like expression.
    """
    from .. import jobparams_naked_end, unix_re
    from ..misc import RelativePath
    from .ordered_dict import OrderedDict

    super(AbstractMassExtract, self).__init__()

    # this fools the derived classes' __setattr__
    self.__dict__.update({'dicttype': dict, '_view': None, 'naked_end': naked_end,
                          'unix_re': unix_re, '_excludes': excludes, 
                          '_cached_extractors': None, 'dynamic': dynamic })
    self.naked_end = jobparams_naked_end if naked_end is None else naked_end
    """ If True and dict to return contains only one item, returns value itself. """
    self.unix_re = unix_re
    """ If True, then all regex matching is done using unix-command-line patterns. """
    self.excludes = excludes
    """ Patterns to exclude. """
    self._cached_extractors = None
    """ List of extration objects. """
    self.dynamic = dynamic
    """ If True chooses a slower but more dynamic caching method. """
    self.dicttype = OrderedDict if ordered else dict
    """ Type of dictionary to use. """
    if path is None: self.__dict__['_rootpath'] = None
    else: self.__dict__['_rootpath']= RelativePath(path, hook=self.uncache)

  @property
  def rootpath(self): 
    """ Root of the directory-tree to trawl for OUTCARs. """
    return self._rootpath.path if self._rootpath is not None else None
  @rootpath.setter
  def rootpath(self, value):
    from ..misc import RelativePath
    if self._rootpath is None:
      self._rootpath = RelativePath(path=value, hook=self.uncache)
    else: self._rootpath.path = value

  def uncache(self): 
    """ Uncache values. """
    self._cached_extractors = None

  @property 
  def excludes(self):
    """ Pattern or List of patterns to ignore. or None.

        :py:attr:`unix_re` determines whether these are unix-command-line like
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

        :param excludes: Pattern or patterns to exclude from output.
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
    if excludes is None or len(excludes) == 0: return self.shallow_copy(excludes=None)
    result = self.shallow_copy(excludes=excludes)
    if self.excludes is not None: result.excludes.extend(self.excludes)
    return result

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

  def __contains__(self, key):
    """ Returns True if key is valid and not empty. """
    from re import compile
    rekey = compile(key)
    for key in self.iterkeys(): 
      if rekey.match(key): return True
    return False

  @staticmethod
  def _translate_regex(pat):
    """ Translates a pattern from unix to re. 

        Compared to fnmatch.translate, doesn't use '.', but rather '[^/]'.
        And doesn't add the tail that fnmatch.translate does.
        Otherwise, code is taked from fnmatch.translate.
    """
    from re import escape
    i, n = 0, len(pat)
    res = ''
    while i < n:
        c = pat[i]
        i = i+1
        if c == '*':
            res = res + '[^/]*'
        elif c == '?':
            res = res + '[^/]'
        elif c == '[':
            j = i
            if j < n and pat[j] == '!':
                j = j+1
            if j < n and pat[j] == ']':
                j = j+1
            while j < n and pat[j] != ']':
                j = j+1
            if j >= n:
                res = res + '\\['
            else:
                stuff = pat[i:j].replace('\\','\\\\')
                i = j+1
                if stuff[0] == '!':
                    stuff = '^' + stuff[1:]
                elif stuff[0] == '^':
                    stuff = '\\' + stuff
                res = '{0}[{0}]'.format(res, stuff)
        else:
            res = res + escape(c)
    return res 

  def _regex_pattern(self, pattern, flags=0):
    """ Returns a regular expression. """
    from re import compile
    from fnmatch import translate
    if self.unix_re: pattern = self._translate_regex(pattern)
    if len(pattern) == 0: return compile("", flags)
    if pattern[-1] in ('/', '\Z', '$'): return compile(pattern, flags)
    return compile(pattern + r"(?=/|\Z)(?ms)", flags)

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
    if self._view is None: return ""
    return self._view
  @view.setter
  def view(self, value): self._view = value

  @property
  def _extractors(self):
    """ Goes through all jobs and collects Extract if available. """
    if self.dynamic:
      if self._cached_extractors is None: self._cached_extractors = self.dicttype()
      result = self.dicttype()
      for name, extract in self.__iter_alljobs__():
        if name not in self._cached_extractors: self._cached_extractors[name] = extract
        result[name] = self._cached_extractors[name]
      return result
    else:
      if self._cached_extractors is not None: return self._cached_extractors
      result = self.dicttype()
      for name, extract in self.__iter_alljobs__(): result[name] = extract
      self._cached_extractors = result
      return result

  def _regex_extractors(self):
    """ Loops through jobs in this view. """
    if self._view is None or self._view == "": 
      for key, value in self._extractors.iteritems(): yield key, value
      return

    regex = self._regex_pattern(self.view)
    if self.excludes is not None: excludes = [self._regex_pattern(u) for u in self.excludes]
    for key, value in self._extractors.iteritems():
      if regex.match(key) is None: continue
      if self.excludes is not None and any(u.match(key) is not None for u in excludes): continue
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
    from .forwarding_dict import ForwardingDict
    assert name in self._attributes, AttributeError("Unknown attribute {0}.".format(name))

    result = self.dicttype()
    for key, value in self.iteritems():
      try: result[key] = getattr(value, name)
      except: result.pop(key, None)
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return ForwardingDict(dictionary=result, naked_end=self.naked_end)

  def __getitem__(self, name):
    """ Returns a view of the current job-dictionary.
    
        .. note:: os.path.normpath_ returns a valid path when descending below
           root, e.g.``normpath('/../../other') == '/other'), so there won't be
           any errors on that account.
    """
    from os.path import normpath, join
    if name[0] != '/': name = join(self.view, name)
    if self.unix_re: name = normpath(name)
    return self.shallow_copy(view=name)

  def __getstate__(self): 
    d = self.__dict__.copy()
    return d

  def __setstate__(self, arg):
    self.__dict__.update(arg)
       
  def shallow_copy(self, **kwargs):
    """ Returns a shallow copy. 
    
        :param kwargs:  Any keyword attribute will modify the corresponding
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
    warn( DeprecationWarning('jobs property is deprecated in favor of keys and iterkeys.'),
          stacklevel=2 )
    return self.keys()

  def iterfiles(self, **kwargs):
    """ Iterates over output/input files. 

        This is rerouted to all extraction objects.
    """
    for job in self.itervalues(): 
      if hasattr(job, 'iterfiles'): 
        for file in job.iterfiles(**kwargs): yield file 

  def __getstate__(self):
    d = self.__dict__.copy()
    if d["_rootpath"] is not None: d["_rootpath"].hook = None
    return d

  def __setstate__(self, arg):
    self.__dict__.update(arg)
    if self._rootpath is not None: self._rootpath.hook = self.uncache
       

class AbstractMassExtractDirectories(AbstractMassExtract):
  """ Propagates extractors from all subdirectories.
  
      Trolls through all subdirectories for calculations with given extraction
      files, and organises results as a dictionary where keys are the name of
      the diretory.

      An class derived from this one should make sure that:
      
      - `Extract` is not none.
      - `__is_calc_dir__` is correctly defined. 
  """
  def __init__(self, path = '.', Extract = None, **kwargs):
    """ Initializes AbstractMassExtractDirectories.
    
    
        :param Extract
            Extraction class to use within each calculation. 
        :param kwargs : dict
            Keyword parameters passed on to :py:meth:`AbstractMassExtract.__init__`.
    """
    # this will throw on unknown kwargs arguments.
    super(AbstractMassExtractDirectories, self).__init__(**kwargs)

    self.__dict__['Extract'] = Extract
    """ Extraction class to use. """

  def __iter_alljobs__(self):
    """ Goes through all directories with a contcar. """
    from os import walk
    from os.path import relpath, join

    for dirpath, dirnames, filenames in walk(self.rootpath, topdown=True, followlinks=True):
      if not self.__is_calc_dir__(dirpath, dirnames, filenames): continue

      result = self.Extract(join(self.rootpath, dirpath))

      result.OUTCAR = self.OUTCAR
      yield join('/', relpath(dirpath, self.rootpath)), result

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
    result = self.__class__(self.rootpath)
    for k, v in self.__dict__.items():
      if k == 'dicttype': result.__dict__[k] = v
      elif k != '_rootpath': result.__dict__[k] = copy(v)
    return result
