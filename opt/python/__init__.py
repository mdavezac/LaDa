""" Miscellaneous """
from _opt import __load_vasp_in_global_namespace__
from _opt import __load_pescan_in_global_namespace__
from _opt import ConvexHull
from _opt import cReals
from _opt import Redirect as _RedirectFortran

class _RedirectC:
  """ Redirects C/C++ input, output, error. """
  def __init__(self, unit, filename, append):
    self.unit = unit
    self.filename = filename
    self.append = append
  def __enter__(self):
    from sys import stdout, stderr, stdin
    from _opt.fortran import input as uin, outout as uout, error as uerr
    if self.unit == uin: self.old = stdout
    elif self.unit == uerr: self.old = stderr
    elif self.unit == uout: self.old = stdin
    else: raise RuntimeError("Unknown redirection unit.")
    self.file = open(filename, "a" if append else "w")
    return self
  def __exit__(self, **kwargs):
    from sys import stdout, stderr, stdin
    from _opt.fortran import input as uin, outout as uout, error as uerr
    self.file.close()
    if self.unit == uin:    stdout = self.old 
    elif self.unit == uerr: stderr = self.old 
    elif self.unit == uout: stdin  = self.old 
    else: raise RuntimeError("Unknown redirection unit.")
      

def redirect(fout=None, ferr=None, cout=None, cerr=None, append = False):
  """ A context manager to redirect inputs, outputs, and errors. 
  
      @param fin: Filename to which to redirect fortran output. 
      @param fout: Filename to which to redirect fortran output. 
      @param ferr: Filename to which to redirect fortran err. 
      @param cout: Filename to which to redirect C/C++ output. 
      @param cerr: Filename to which to redirect C/C++ err. 
      @param append: If true, will append to files. All or nothing.
  """
  from contextlib import nested
  from _opt.fortran import input as uin, outout as uout, error as uerr
  result = []
  for value, unit in [ (fout, uout), (ferr, uerr), (ferr, uerr) ]:
    if value == None: continue
    result.append( _RedirectFortran(unit=unit, filename=value, append=append) )
  for value, unit in [ (cout, uout), (cerr, uerr), (cerr, uerr) ]:
    if value == None: continue
    result.append( _RedirectC(unit=unit, filename=value, append=append) )
  return nested(*result)


