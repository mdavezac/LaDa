""" Miscellaneous """
import _opt 
__load_vasp_in_global_namespace__ = _opt.__load_vasp_in_global_namespace__
__load_pescan_in_global_namespace__ = _opt.__load_pescan_in_global_namespace__
cReals = _opt.cReals
_RedirectFortran = _opt._RedirectFortran

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
    self.file = open(filename if len(filename) else "/dev/null", "a" if append else "w")
    return self
  def __exit__(self, **kwargs):
    from sys import stdout, stderr, stdin
    from _opt.fortran import input as uin, outout as uout, error as uerr
    self.file.close()
    if self.unit == uin:    stdout = self.old 
    elif self.unit == uerr: stderr = self.old 
    elif self.unit == uout: stdin  = self.old 
    else: raise RuntimeError("Unknown redirection unit.")
      

def redirect(fout=None, ferr=None, fin=None, cout=None, cerr=None, cin=None, append = False):
  """ A context manager to redirect inputs, outputs, and errors. 
  
      @param fout: Filename to which to redirect fortran output. 
      @param ferr: Filename to which to redirect fortran err. 
      @param fin: Filename to which to redirect fortran input. 
      @param cout: Filename to which to redirect C/C++ output. 
      @param cerr: Filename to which to redirect C/C++ err. 
      @param cin: Filename to which to redirect C/C++ input. 
      @param append: If true, will append to files. All or nothing.
  """
  from contextlib import nested
  from _opt.fortran import input as uin, outout as uout, error as uerr
  result = []
  for value, unit in [ (fout, uout), (ferr, uerr), (fin, uin) ]:
    if value == None: continue
    result.append( _RedirectFortran(unit=unit, filename=value, append=append) )
  for value, unit in [ (cout, uout), (cerr, uerr), (cin, uin) ]:
    if value == None: continue
    result.append( _RedirectC(unit=unit, filename=value, append=append) )
  return nested(*result)


