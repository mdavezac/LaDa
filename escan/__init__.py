""" Interface module for ESCAN. """
__docformat__ = "restructuredtext en"
__all__ = [ "Extract", 'MassExtract', "bandgap", "extract_bg",
            "dipole_matrix_element", "band_structure", "call_escan",
            "Escan", "folded_spectrum", "all_electron", "soH", 
            "nonlocalH", "localH", "AtomicPotential", "band_structure",
            "extract_bg", 'ExtractBS', 'KEscan', 'KPoints', 'KGrid', 
            'ReducedKGrid', 'ReducedKDensity', 'BPoints', 'ReducedBPoints', 
            'soH', 'nonlocalH', 'localH', 'folded_spectrum', 'all_electron',
            'read_input', 'exec_input']

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
from _bandstructure import band_structure, Extract as ExtractBS, BPoints, ReducedBPoints
from _bandgap import bandgap, extract as extract_bg
from _extract import Extract, MassExtract
from _extract import Extract as _EscanExtract
from functional import Escan, folded_spectrum, all_electron
from _potential import soH, nonlocalH, localH, AtomicPotential
from kescan import KEscan, Extract as KExtract
from kpoints import KGrid, ReducedKGrid, KPoints, ReducedKDensity


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


def exec_input(filepath = "input.py", namespace = None):
  """ Executes an input script including namespace for escan/vff. """ 
  from ..vff import Vff
  from ..opt import exec_input
  from . import Escan, soH, nonlocalH, localH, folded_spectrum, all_electron

  dictionary = { "Vff": Vff, "Escan": Escan, "soH": soH, \
                 "nonlocalH": nonlocalH, "localH": localH, \
                 "folded_spectrum": folded_spectrum, "all_electron": all_electron}
  if namespace != None: dictionary.update(namespace)
  return exec_input(filepath, dictionary)

def read_input(filepath = "input.py", namespace = None, name=None):
  """ Reads an input file including namespace for escan/vff. """ 
  from ..vff import Vff
  from ..opt import read_input
  from . import Escan, soH, nonlocalH, localH, folded_spectrum, all_electron

  dictionary = { "Vff": Vff, "Escan": Escan, "soH": soH, \
                 "nonlocalH": nonlocalH, "localH": localH, \
                 "folded_spectrum": folded_spectrum, "all_electron": all_electron}
  if namespace != None:
    dictionary.update(namespace)
    if name == None and hasattr(namespace, '__name__'): namee = namespace.__name__
  return read_input(filepath, dictionary)
