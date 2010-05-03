""" Miscellaneous """
import _opt 
__load_vasp_in_global_namespace__ = _opt.__load_vasp_in_global_namespace__
__load_pescan_in_global_namespace__ = _opt.__load_pescan_in_global_namespace__
cReals = _opt.cReals
_RedirectFortran = _opt._RedirectFortran
streams = _opt._RedirectFortran.fortran

class _RedirectC:
  """ Redirects C/C++ input, output, error. """
  def __init__(self, unit, filename, append):
    self.unit = unit
    self.filename = filename
    self.append = append
  def __enter__(self):
    import sys
    if self.unit == streams.input:    self.old = sys.stdin
    elif self.unit == streams.error:  self.old = sys.stderr
    elif self.unit == streams.output: self.old = sys.stdout
    else: raise RuntimeError("Unknown redirection unit.")
    self.file = open(self.filename if len(self.filename) else "/dev/null", "a" if self.append else "w")
    if self.unit == streams.input:    sys.stdin  = self.file
    elif self.unit == streams.error:  sys.stderr = self.file
    elif self.unit == streams.output: sys.stdout = self.file
    else: raise RuntimeError("Unknown redirection unit.")
    return self
  def __exit__(self, *wargs):
    import sys 
    if self.unit == streams.input:    sys.stdin  = self.old 
    elif self.unit == streams.error:  sys.stderr = self.old 
    elif self.unit == streams.output: sys.stdout = self.old 
    else: raise RuntimeError("Unknown redirection unit.")
    self.file.close()
      

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
  result = []
  for value, unit in [ (fout, streams.output), (ferr, streams.error), (fin, streams.input) ]:
    if value == None: continue
    result.append( _RedirectFortran(unit=unit, filename=value, append=append) )
  for value, unit in [ (cout, streams.output), (cerr, streams.error), (cin, streams.input) ]:
    if value == None: continue
    result.append( _RedirectC(unit=unit, filename=value, append=append) )
  return nested(*result)


