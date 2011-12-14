""" Interface module for ESCAN. """
__docformat__ = "restructuredtext en"
__all__ = [ "Extract", 'MassExtract', "bandgap", "extract_bg", 'ldos',
            "Escan", "folded_spectrum", "all_electron", "soH", 'BandGap',
            "nonlocalH", "localH", "AtomicPotential", "extract_bg", 'KEscan', 'KPoints', 
            'KGrid', 'ReducedKGrid', 'ReducedKDensity', 'soH', 'nonlocalH', 'localH', 
            'folded_spectrum', 'all_electron', 'read_input', 'exec_input', 'KExtract',  
            'majority_representation', 'BPoints', 'ReducedBPoints', 'plot_bands', 'plot_alloybands',
            'fftmesh', 'EMassFunctional', 'EMassExtract', 'InnerBPoints', 'ReducedInnerBPoints' ]

from _bandstructure import plot_bands, BPoints, ReducedBPoints,\
                           plot_alloybands, InnerBPoints, ReducedInnerBPoints
from _bandgap import bandgap, extract as extract_bg, Functional as BandGap
from emass import Functional as EMassFunctional, Extract as EMassExtract
from _extract import Extract
from _massextract import MassExtract
from _potential import soH, nonlocalH, localH, AtomicPotential
from _methods import majority_representation
from functional import Escan, folded_spectrum, all_electron
from kescan import KEscan, Extract as KExtract
from kpoints import KGrid, ReducedKGrid, KPoints, ReducedKDensity
import ldos
import fftmesh


def exec_input(script, namespace = None):
  """ Executes an input script including namespace for escan/vff. """ 
  from ..opt import exec_input as opt_exec_input
  from .. import vff

  dictionary = {}
  for key in vff.__all__: dictionary[key] = getattr(vff, key)
  for key in __all__: dictionary[key] = globals()[key]
  if namespace is not None: dictionary.update(namespace)
  return opt_exec_input(script, dictionary)

def read_input(filepath = "input.py", namespace = None):
  """ Reads an input file including namespace for escan/vff. """ 
  from ..opt import read_input as opt_read_input
  from .. import vff

  dictionary = {}
  for key in vff.__all__: dictionary[key] = getattr(vff, key)
  for key in __all__: dictionary[key] = globals()[key]
  if namespace is not None: dictionary.update(namespace)
  return opt_read_input(filepath, dictionary)

def effective_mass_tensor(escan, structure, outdir=None, comm=None, attenuate=None, degeneracy=None, **kwargs):
  """ Computes effective mass tensor using dipoles and f-sum rule. 
  
      :Parameters:
        escan : escan.Escan
          Straight off escan functional. Calculations must be full diagonalization.
        structure : lada.crystal.Structure
          Crystal structure for which to perform calculations.
        outdir : str or None
          Output directory. Defaults to current directory.
        comm : lada.mpi.communicator.
          MPI communicator with which to perform calculation. 

      :return: a numpy array with units of an electronic mass, with the 
               dimensions of the input directions, except for the last axis,
               which is of the same length as the number of computed
               eigenvalues.
  """
  if kwargs.get('eref', escan.eref) != None: 
    raise ValueError("f-sum rule computation of effective mass only done with direct diagonalization.")

  # compute escan wavefunctions.
  result = escan(structure, outdir, comm, **kwargs)

  # compute effective mass tensor.
  ekwargs = {}
  if attenuate != None: ekwargs['attenuate'] = attenuate
  if degeneracy != None: ekwargs['degeneracy'] = degeneracy
  return result.effective_mass_tensor(**ekwargs)

effective_mass_tensor.Extract = Extract
""" Extraction object for f-sum rule effective masses. """
