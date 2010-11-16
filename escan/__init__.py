""" Interface module for ESCAN. """
__docformat__ = "restructuredtext en"
__all__ = [ "Extract", 'MassExtract', "bandgap", "extract_bg",\
            "dipole_matrix_element", "band_structure", "call_escan",\
            "Escan", "folded_spectrum", "all_electron", "soH", \
            "nonlocalH", "localH", "AtomicPotential", "band_structure",\
            "extract_bg", 'ExtractBS', 'KEscan', 'KPoints', 'KGrid', \
            'ReducedKGrid', 'ReducedKDensity' ]

from ..opt import __load_escan_in_global_namespace__
from lada import lada_with_mpi
if lada_with_mpi:
  if __load_escan_in_global_namespace__:
    from DLFCN import RTLD_NOW as _RTLD_NOW, RTLD_GLOBAL as _RTLD_GLOBAL
    from sys import getdlopenflags as _getdlopenflags, setdlopenflags as _setdlopenflags
    flags = _getdlopenflags()
    _setdlopenflags(_RTLD_NOW|_RTLD_GLOBAL)
    import _escan
    _setdlopenflags(flags)
  else: import _escan
from ..opt.decorators import broadcast_result
from _bandstructure import band_structure, Extract as ExtractBS
from _bandgap import bandgap, extract as extract_bg
from _extract import Extract, MassExtract
from _extract import Extract as _EscanExtract
from functional import Functional as Escan
from kescan import KEscan, KGrid, ReducedKGrid, KPoints, ReducedKDensity


def call_escan(comm, atom="atom_input", pot="pot_input", escan="escan_input"):
  """ Calls escan functional in current directory.

      :Parameters:
        comm : boost.mpi.communicator
          Processes on which to execute.
        atom
          file with atomic input (from vff). 
        pot
          file with input to potential generation.
        escan
          file with input to escan itself.

      Before calling the functional, the files are propagated such that each
      proc can read its own. What is more, the inputs are modified such that
      L{pot} does refer to L{atom}, and L{escan} to both atomic config and
      potential files.
  """
  from os import remove
  if lada_with_mpi:
    from boost.mpi import world, broadcast
    from _escan import _call_escan, _call_genpot
  else:
    class world:
      rank = 0
      def barrier(self): pass
    comm = world
    def broadcast(comm, value, *args, **kwargs): return value
  
  # escan is free for all, with every proc reading at the same time.
  # hence it must be from different files.
  atominput = "%s.%i" % (atom, world.rank)
  potinput = "%s.%i" % ("pot.input", world.rank)
  potoutput = "%s.%i" % ("pot.output", world.rank)
  escaninput = "%s.%i" % ("escan_input", world.rank)

  # propagates atomic input first
  buffer = None
  if comm.rank == 0:
    with open(atom, "r") as file: buffer = file.read()
  buffer = broadcast(comm, buffer, 0)
  with open(atominput, "w") as file: file.write(buffer)

  # propagates pot input second
  buffer = None
  if comm.rank == 0:
    with open(pot, "r") as file: buffer = file.readlines()
    buffer[0] = buffer[0].replace( buffer[0].split()[0], atom )
  buffer = broadcast(comm, buffer, 0)
  with open(potinput, "w") as file: file.writelines(buffer)

  # propagates pot input second
  buffer = None
  if comm.rank == 0:
    with open(escan, "r") as file: buffer = file.readlines()
    buffer[0] = buffer[0].replace( buffer[0].split()[1], potoutput )
    buffer[12] = buffer[12].replace( buffer[12].split()[1], atominput )
  buffer = broadcast(comm, buffer, 0)
  with open(escaninput, "w") as file: file.writelines(buffer)

  #  calls escan at last. 
  if lada_with_mpi: 
    _call_genpot(comm)
    _call_escan(comm)
  else: 
    raise RuntimeError('Cannot call escan if lada_with_mpi == False.')

  comm.barrier()
  remove(atominput)
  remove(potinput)
  remove(escaninput)

folded_spectrum = 0
""" Folded spectrum method. """
all_electron = 1
""" All electron method. """
localH = 0
""" Local hamiltonian. """
nonlocalH = 1
""" Non-local hamiltonian. """
soH = 2
""" Spin-orbit hamiltonian. """

class AtomicPotential(object):
  """ Holds parameters to atomic potentials. """
  def __init__(self, path, nonlocal=None, s=None, p=None, d=None, pnl=None, dnl=None):
    from ..opt import RelativeDirectory

    self._filepath = RelativeDirectory(path=path)
    """ Private path to pseudopotential file. 
    
        Path is a relative directory for added transferability from computer to
        computer.
    """
    self._nonlocal = None if nonlocal == None else RelativeDirectory(nonlocal)
    """ Private path to non-local part, or None. 
    
        Path is a relative directory for added transferability from computer to
        computer.
    """
    self.s =  s if s != None else 0
    """ s parameter """
    self.p =  p if p != None else 0
    """ p parameter """
    self.d =  d if d != None else 0
    """ d parameter """
    self.pnl = pnl if pnl != None else 0
    """ p non-local parameter """
    self.dnl = dnl if dnl != None else 0
    """ d non-local parameter """

  @property
  def filepath(self):
    """ Path to pseudopotential file. """
    return self._filepath.path
  @filepath.setter
  def filepath(self, value): self._filepath.path = value

  @property
  def nonlocal(self):
    """ Path to pseudopotential file. """
    return self._nonlocal.path
  @nonlocal.setter
  def nonlocal(self, value): self._nonlocal.path = value

  def __repr__(self):
    """ Tuple representation of self. """
    from os.path import relpath
    result = '"{0}"'.format(self._filepath.unexpanded)
    if self.nonlocal == None: result += ", None"
    else: result += ', "%s"' % (self._nonlocal.unexpanded)
    result += ", %f, %f, %f, %f, %f" % (self.s, self.p, self.d, self.pnl, self.dnl)
    return result

  @broadcast_result(key=True)
  def get_izz(self, comm = None):
    """ Returns izz string greped from pseudopotential file. """
    with open(self.filepath, "r") as file:
      return file.readline().split()[1]


def read_input(filepath = "input.py", namespace = None):
  """ Reads an input file including namespace for escan/vff. """ 
  from ..jobs import JobDict
  from ..vff import Vff
  from ..opt import read_input
  from . import Escan, soH, nonlocalH, localH, folded_spectrum, all_electron

  dictionary = { "Vff": Vff, "Escan": Escan, "soH": soH, \
                 "nonlocalH": nonlocalH, "localH": localH, \
                 "folded_spectrum": folded_spectrum, "all_electron": all_electron}
  if namespace != None: dictionary.update(namespace)
  return read_input(filepath, dictionary)
