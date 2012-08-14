__docformat__ = "restructuredtext en"
__all__ = ['Transforms', 'supercells']
from .transforms import Transforms

def supercells(lattice, sizerange):
  """ Generator over supercells for given size range.
  
      :params lattice: 
         Back-bone lattice
      :type lattice: py:attr:`~lada.crystal.Structure`
      :param sizerange: 
         List of sizes for which to perform calculations, in number of
         unit-cells per supercell.
      :type sizerange: integer sequence
      :returns:
          dictionary where each key is a size and each element a list of
          inequivalent supercells of that size
  """
  from itertools import product
  from numpy import dot, mod, equal, array, all
  from numpy.linalg import inv
  from ..crystal import space_group
  from .cppwrappers import is_integer

  sizerange = sorted([k for k in sizerange if k > 0])
  results = {}
  for n in sizerange: results[n] = []
  mink = min(sizerange)
  maxk = max(sizerange)
  cell = lattice.cell
  invcell = inv(cell)
  spacegroup = space_group(lattice)
  for i in xrange(len(spacegroup)):
    spacegroup[i] = dot(invcell, dot(spacegroup[i][:3], cell))

  def isthere(sc, l):
    """ Check if a known supercell """
    isc = inv(sc)
    for op in spacegroup:
      op2 = dot(isc, op)
      if any( is_integer(dot(op2, u)) for u in l): return True
    return False

  supercell = array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype='float64')
  for a in range(1, maxk+1):
    maxb = maxk // a 
    if maxb % a != 0: maxb += 1 
    supercell[0,0] = a
    for b in xrange(1, maxb+1): 
      maxc = maxk//(a*b)
      if maxk % (a*b) != 0: maxc += 1
      supercell[1,1] = b
      for c in xrange(max(mink//(a*b), 1), maxc+1):
        n = a * b * c
        if n not in sizerange: continue
        supercell[2, 2] = c
        for d, e, f in product(xrange(b), xrange(c), xrange(c)):
          supercell[1, 0] = d
          supercell[2, 0] = e
          supercell[2, 1] = f
          if not isthere(supercell, results[n]): 
            results[n].append(supercell.copy())

  return results




def hf_groups(lattice, sizerange):
  """ Generator over supercells for given size range.
  
      :params lattice: 
         Back-bone lattice
      :type lattice: py:attr:`~lada.crystal.Structure`
      :param sizerange: 
         List of sizes for which to perform calculations, in number of
         unit-cells per supercell.
      :type sizerange: integer sequence
      :yields: two tuple where the first 
  """
  pass
