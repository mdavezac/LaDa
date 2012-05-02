""" Holds exceptions declared by LaDa. """

class root(Exception):
  """ Root for all LaDa exceptions. """
  pass

class input(root):
  """ Root for all input LaDa exceptions. """
  pass

class out_of_range(root):
  """ Root for all out-of-range LaDa exceptions. """
  pass

class internal(root, RuntimeError):
  """ Root for all internal (cpp) LaDa exceptions. """
  pass

class infinite_loop(root):
  """ Root for all infinite-loops LaDa exceptions. """
  pass

class ValueError(root, ValueError):
  """ Root for all ValueError LaDa exceptions. """
  pass

class KeyError(root, KeyError):
  """ Root for all KeyError LaDa exceptions. """
  pass

class AttributeError(root, AttributeError):
  """ Root for all AttributeError LaDa exceptions. """
  pass

class IndexError(root, IndexError):
  """ Root for all IndexError LaDa exceptions. """
  pass

class TypeError(root, TypeError):
  """ Root for all TypeError LaDa exceptions. """
  pass

class NotImplementedError(root, NotImplementedError):
  """ Root for all NotImplementedError LaDa exceptions. """
  pass

class ImportError(root, ImportError):
  """ Root for all ImportError LaDa exceptions. """
  pass

class IOError(root, IOError):
  """ Root for all ImportError LaDa exceptions. """
  pass

class Math(root):
  """ Root of math exceptions. """
  pass
class singular_matrix(Math):
  """ Singular matrix. """
  pass
