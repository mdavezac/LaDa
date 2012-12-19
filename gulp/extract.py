""" Sub-package for gulp extraction object. """
__docformat__ = "restructuredtext en"
__all__ = ['Extract']
from ..tools import make_cached
from ..tools.extract import search_factory, AbstractExtractBase
OutputSearchMixin = search_factory('OutputSearchMixin', 'stdout', __name__, 'gulp.gout')

class ExtractBase(object):
  """ Greps some stuff from GULP output file. """
  def __init__(self):
    """ Initializes extraction class. """
    super(ExtractBase, self).__init__()

  @property
  @make_cached
  def optimize(self):
    """ True if optimizing structure. """
    return self._find_first_STDOUT('\*\s*optimise\s*\-\s*perform') is not None
  @property
  @make_cached
  def conp(self):
    """ True if optimizing at constant pressure. """
    return self._find_first_STDOUT('\*\s*conp\s*\-\s*constant') is not None
  @property
  @make_cached
  def qeq(self):
    """ True if using QEe. """
    return self._find_first_STDOUT('\*\s*QEq\s*\-\s*electroneg') is not None
  @property
  @make_cached
  def cell(self):
    """ Greps output cell. """
    # static case
    if not self.optimize: return self.input_cell

    from re import M
    from numpy import array
    from quantities import angstrom
    from ..error import GrepError
    expr = '\s*Final\s+Cartesian\s+lattice\s+vectors\s*\(Angstroms\)\s*:\s*'   \
           '\\n\s*\\n'                                                         \
           '(?:\s*(?:(\S*)\s*){3}\s*\\n){3}'
    result = self._find_first_STDOUT(expr, M)
    if result is None: raise GrepError('Could not grep cell')
    return array(result.group(0).split()[-9:], dtype='float64')                \
                .reshape(3,3) * angstrom
  
  @property
  @make_cached
  def input_cell(self):
    """ Greps input cell. """
    from re import M
    from numpy import array
    from quantities import angstrom
    from ..error import GrepError
    expr = '\s*Cartesian\s+lattice\s+vectors\s*\(Angstroms\)\s*:\s*'           \
           '\\n\s*\\n'                                                         \
           '(?:\s*(?:(\S*)\s*){3}\s*\\n){3}'
    result = self._find_first_STDOUT(expr, M)
    if result is None: raise GrepError('Could not grep cell')
    return array(result.group(0).split()[-9:], dtype='float64')                \
                .reshape(3,3) * angstrom

  @property 
  @make_cached
  def structure(self):
    """ Greps output structure. """
    # static case
    if not self.optimize: return self.input_structure

    from re import compile
    from itertools import chain
    from quantities import e, angstrom
    from ..crystal import Structure
    from ..error import GrepError
    regex = compile('Final asymmetric unit coordinates') 
    with self.__stdout__() as file:
      # find beginning of final structure
      found = None
      for line in file:
        found = regex.search(line)
        if found is not None: break
      if found is None: raise GrepError('Could not grep atoms.')
      # find first atom
      regex = compile('(\d+)\s+([A-Z][a-z]?(?:\d+)?)\s+\S(\s+\S+){5}')
      for line in file:
        found = regex.search(line)
        if found: break
      if found is None: raise GrepError('Could not grep atoms.')
      # create result
      result = Structure()
      result.cell = self.cell
      for line in chain([line], file):
        found = regex.search(line)
        if found is None: break
        data = line.split()
        result.add_atom( pos=[float(i) for i in data[3:6]],
                         type=data[1], label=int(data[0]) )
        radius = float(data[6])
        if abs(radius) > 1e-8: result[-1].radius = radius * angstrom
        charge = float(data[7])
        if abs(charge) > 1e-8: result[-1].charge = charge * e
      # add more info
      try: result.energy = self.energy
      except: pass
      return result

  @property 
  @make_cached
  def input_structure(self):
    """ Greps input structure. """
    from re import compile
    from itertools import chain
    from quantities import e
    from ..crystal import Structure
    from ..error import GrepError
    regex = compile('Fractional coordinates of asymmetric unit :')
    with self.__stdout__() as file:
      # find beginning of final structure
      found = None
      for line in file:
        found = regex.search(line)
        if found is not None: break
      if found is None: raise GrepError('Could not grep atoms.')
      # find first atom
      regex = compile('(\d+)\s+([A-Z][a-z]?(?:\d+)?)\s+\S(\s+\S+){5}')
      for line in file:
        found = regex.search(line)
        if found: break
      if found is None: raise GrepError('Could not grep atoms.')
      # create result
      result = Structure()
      result.cell = self.input_cell
      for line in chain([line], file):
        found = regex.search(line)
        if found is None: break
        data = line.replace('*', '').split()
        result.add_atom( pos=[float(i) for i in data[3:6]],
                         type=data[1], label=int(data[0]) )
        charge = float(data[6])
        if abs(charge) > 1e-8: result[-1].charge = charge * e
        occupancy = float(data[7])
        if abs(occupancy) > 1e-8: result[-1].occupancy = occupancy
      return result

  @property
  @make_cached
  def energy(self):
    """ Greps energy from result. """
    from re import M
    from quantities import eV
    from ..error import GrepError
    regex = 'Total\s+lattice\s+energy\s*:\s*\\n'                               \
            '\s*Primitive\s*unit\s*cell\s*=\s*(\S+)\s*eV'
    result = self._find_first_STDOUT(regex, M)
    if result is None: raise GrepError('Could not find energy.')
    return float(result.group(1)) * eV

  @property
  def _finished(self):
    """ True if job finished. """
    return self._find_first_STDOUT('Job Finished at') is not None



class Extract(AbstractExtractBase, OutputSearchMixin, ExtractBase):
  """ Extracts GULP data from an output file. """
  def __init__(self, directory=None, **kwargs):
    """ Initializes extraction object. 
    
        :param directory: 
          Directory where the OUTCAR resides. 
          It may also be the path to an OUTCAR itself, if the file is not
          actually called OUTCAR.
    """
    from os.path import exists, isdir, basename, dirname
    from ..misc import RelativePath
       
    self.STDOUT = 'gulp.gout'
    """ Name of file to grep. """
    if directory is not None:
      directory = RelativePath(directory).path
      if exists(directory) and not isdir(directory):
        self.STDOUT = basename(directory)
        directory = dirname(directory)
    AbstractExtractBase.__init__(self, directory)
    ExtractBase.__init__(self)
    OutputSearchMixin.__init__(self)
  @property
  def success(self): 
    """ True if calculations ran to completion. """
    try: return self._finished
    except: return False
