""" Contains basic data type and methods for crystal structure and lattices.
    The basic types are imported in this namespace from (and described in)
    module L{_crystal}.
"""
from _crystal import *

def deform_kpoint(kpoint, ideal, relaxed):
  """ Deform kpoints from ideal cell to relaxed cell. 

      @param kpoint: The kpoint to deform in cartesian coordinates.
      @type kpoint: numpy array
      @param ideal: The original (real-space) cell, as an ideal supercell of the lattice.
      @type ideal: numpy 2d-array
      @param relaxed: The relaxed (real-space) cell.
      @type relaxed: numpy 2d-array
      @return: the kpoint deformed from the ideal reciprocal cell to the
               relaxed reciprocal cell, in cartesian coordinates.
  """
  from numpy import dot, matrix
  from numpy.linalg import inv
  k = dot(ideal.T, kpoint)
  for i in range(3):
    k[i] = float(k[i]) - float( int(k[i]) )
    if k[i] < 0e0: k[i] += 1e0
    if k[i] > 1e0-1e-6: k[i] = 0e0
  return dot(inv(relaxed.T), k.T)
  

