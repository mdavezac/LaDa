""" Mass extraction object for ladabase. """
__all__ = ['MassExtract']
from .extract import VaspExtract

class MassExtract(object): 
  """ Propagates extraction methods from different jobs. """
  Extract = staticmethod(VaspExtract)
  """ Extraction class. """
  def __init__(self, filters=None, namespace=None, **kwargs):
    """ Initializes extraction object. 

        :Parameters:
          filters : list of str or None
            List of filters through which to view the database.
          namespace : dict/namespace
            Names which may be present when evaluating the filters.

        :Kwarg naked_end: True if should return value rather than dict when only one item.
    """
    from .. import naked_end
    object.__init__(self)

    self.naked_end = kwargs.pop('naked_end', naked_end)
    """ If True and dict to return contains only one item, returns value itself. """
    self.filters = [] if filters == None else filters
    """ The filters which jobs should match. """
    self._cached_extractors = None
    """ List of extration objects. """
    self.namespace = {} if namespace == None else namespace
    """ Namespace with which to perform evaluation. """

  def uncache(self): 
    """ Uncache values. """
    self._cached_extractors = None

  def iteritems(self):
    """ Iterates through all extraction objects and names. """
    for name, job in self._filtered_extractors(): yield name, job
  def items(self):
    """ Iterates through all extraction objects and names. """
    return [(name, job) for name, job in self.iteritems()]
    
  def itervalues(self):
    """ Iterates through all extraction objects. """
    for name, job in self._filtered_extractors(): yield job
  def values(self):
    """ Iterates through all extraction objects. """
    return [job for job in self.itervalues()]

  def iterkeys(self):
    """ Iterates through all extraction objects. """
    for name, job in self._filtered_extractors(): yield name
  def keys(self):
    """ Iterates through all extraction objects. """
    return [name for name in self.iterkeys()]
  
  def __iter__(self):
    """ Iterates through all job names. """
    for name, job in self.iteritems(): yield name
  def __len__(self): 
    """ Returns length of output dictionary. """
    return len(self.keys())

  def __contains__(self, filter):
    """ Returns True if filter is valid and not empty. """
    for extract in self.itervalues(): 
      if self._eval_filter(filter, extract): return True
    return False

  def _eval_filter(self, filter, extract):
    """ Returns a regular expression. """
    from math import pi 
    from os import environ
    from numpy import array, matrix, dot, sqrt, abs, ceil
    from numpy.linalg import norm, det
    from .. import physics, crystal

    namespace = { "environ": environ, "pi": pi, "array": array, "matrix": matrix, "dot": dot,
                  "norm": norm, "sqrt": sqrt, "ceil": ceil, "abs": abs,  "det": det,
                  "physics": physics, 'crystal': crystal, 'id': str(extract._id) } 
    namespace.update(self.namespace)
    namespace["extract"] = extract

    return eval(filter, namespace)

  def __iter_alljobs__(self):
    """ Goes through database. """
    from .misc import get_ladabase
    for doc in get_ladabase().files.find():
      yield str(doc['_id']), self.Extract(doc)

  @property
  def _extractors(self):
    """ Goes through all jobs and collects Extract if available. """
    if self._cached_extractors != None: return self._cached_extractors
    result = {}
    for name, extract in self.__iter_alljobs__(): result[name] = extract
    self._cached_extractors = result
    return result

  def _filtered_extractors(self):
    """ Loops through jobs in this view. """
    for key, value in self._extractors.iteritems():
      dothis = True
      for filter in self.filters:
        if not self._eval_filter(filter, value): dothis = False; break
      if dothis: yield key, value

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
    return list(set(results)) if len(self) != 1 else list(set(results) | set('outcar'))

  def __getattr__(self, name): 
    """ Returns extracted values. """
    from ..jobs.forwarding_dict import ForwardingDict

    assert name in self._attributes, AttributeError("Unknown attribute {0}.".format(name))

    result = {}
    for key, value in self.iteritems():
      try: result[key] = getattr(value, name)
      except: result.pop(key, None)
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return ForwardingDict(result, naked_end=self.naked_end)

  def __getitem__(self, filters):
    """ Returns a view of the current job-dictionary. """
    if isinstance(filters, str): filters = [filters]
    result = self.__class__(namespace = self.namespace)
    extractor = self.__class__(namespace = self.namespace, filters=filters)
    extractor._cached_extractors = self._cached_extractors
    result._cached_extractors = {}
    for key, value in extractor.iteritems(): result._cached_extractors[key] = value
    return result
