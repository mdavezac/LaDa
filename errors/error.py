""" Holds exceptions declared by Pylada. """

class root(Exception):
  """ Root for all Pylada exceptions. """
  pass

class input(root):
  """ Root for all input Pylada exceptions. """
  pass

class out_of_range(root):
  """ Root for all out-of-range Pylada exceptions. """
  pass

class internal(root, RuntimeError):
  """ Root for all internal (cpp) Pylada exceptions. """
  pass

class infinite_loop(root):
  """ Root for all infinite-loops Pylada exceptions. """
  pass

class ValueError(root, ValueError):
  """ Root for all ValueError Pylada exceptions. """
  pass

class KeyError(root, KeyError):
  """ Root for all KeyError Pylada exceptions. """
  pass

class AttributeError(root, AttributeError):
  """ Root for all AttributeError Pylada exceptions. """
  pass

class IndexError(root, IndexError):
  """ Root for all IndexError Pylada exceptions. """
  pass

class TypeError(root, TypeError):
  """ Root for all TypeError Pylada exceptions. """
  pass

class NotImplementedError(root, NotImplementedError):
  """ Root for all NotImplementedError Pylada exceptions. """
  pass

class ImportError(root, ImportError):
  """ Root for all ImportError Pylada exceptions. """
  pass

class IOError(root, IOError):
  """ Root for all ImportError Pylada exceptions. """
  pass

class Math(root):
  """ Root of math exceptions. """
  pass
class singular_matrix(Math):
  """ Singular matrix. """
  pass
class interactive(input):
  """ Interactive usage error. """
  pass
class GrepError(AttributeError):
  """ Raised when property could not be grepped from some OUTCAR. """
  pass
class ConfigError(input):
  """ Some sort of Pylada configuration error. """

class ExternalRunFailed(root):
  """ Thrown when an external run has failed. """
