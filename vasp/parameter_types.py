""" Paremeters types for use as attributes in Incar """
class Standard(object):
  """ A standard parameter in the form of key/value pair. """
  def __init__(self, key, value, validity=None, doc=None):

    self.key = key
    self.value = value
    if doc != None: self.__doc__ = doc

  def __set__(self, parent, value):
    """ Sets value of key/value pair """
    if self.validity != None:
      assert self.validity(value), \
             "Value %s = %s is invalid.\n%s\n"  % (self.key, value, self.__doc__)
    self.value = value

  def __get__(self, parent, parenttype=None):
    """ Gets value of key/value pair """
    return self.value

  def __str__(self):
    """ Prints in INCAR format. """
    return "%s = %s" % (self.key, self.value)

  def add_to_doc(self, parent):

    lock = "_LockDocOfObject%s" % (self.__name__)
    if hasattr(parent.__class__, lock): return 
    setattr(parent, lock, None)
    parent.__doc__ += "self.%s:\n%s" % (self.__name__, self.__doc__)

class NoPrintStandard(Standard):
  """ Does not print if value is the default given on initialization. """
  def __init__(self, key, value, validity=None, doc=None):
    Standard.__init__(self, key, value, validity, doc)
    self.default = value
  def __str__(self):
    if self.default == self.value: 
      return "# %s = VASP default." % (self.key)
    return Standard.__str__(self)

