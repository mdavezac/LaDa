""" Module containing different k-mesh schemes """
class Gamma(object):
  """ K-mesh with a single point: Gamma. """
  def __call__(self, vasp):
    """ Returns a string which is the text of the KPOINTS file """
    return """Gamma only\n0\nGamma\n1 1 1\n0 0 0\n"""

class Density(object):
  """ Contains vasp kpoints parameters. """
  def __init__(self, offset=(0.5, 0.5, 0.5), length=70):
    """ Initializes density.
          - offset is with respect to Gamma, eg (0, 0, 0) means gamma centered.
          - length is the target number of points per 2\pi/Angstrom.
    """

    self.offset = offset
    self.length = length

  def __call__(self, vasp):
    """ Returns a string which is the text of the KPOINTS file """
    from lada.crystal.gruber import Reduction
    from lada.crystal import get_point_group_symmetries as point_group
    import numpy as np
    from math import fabs, floor
    result = "Automatic k-mesh generation\n0\nGamma\n"

    reduction = Reduction()
    cell = np.matrix( [[ vasp.system.cell[j,i] for j in range(3)] for i in range(3)],
                      dtype="float64" )
    assert fabs(np.linalg.det(cell)) > 1e-12, ValueError
    # gets reciprocal cell
    cell = cell.I
    # gets Niggli reciprocal cell
    cell = reduction(cell)
    # gets original point-group symmetries
    pointgroup = point_group(cell)

    # creates trial-grid
    n = [[0, 1+int(floor(max([1e0, np.linalg.norm(cell[:,0]) * float(self.length)])))],
         [1, 1+int(floor(max([1e0, np.linalg.norm(cell[:,1]) * float(self.length)])))], 
         [2, 1+int(floor(max([1e0, np.linalg.norm(cell[:,2]) * float(self.length)])))] ]

    grid  = cell.copy()
    for i in range(3): grid[:,n[i][0]] *= 1e0 / float(n[i][1] )
    if vasp.isym == "off" or Density._goodsyms(grid, pointgroup): # if symmetry are preserved. 
      for i in range(3): result += "%18.12f %18.12f %18.12f\n" % tuple(grid[:,i])
      result += "  %18.12f %18.12f %18.12f\n" % (self.offset[0], self.offset[1], self.offset[2])
      return result

    # if symmetries are not preserved, then figure out smallest diagonal cell
    # I'm not sure this work at all. One solution is to use n[0][1] == n[1][1] == n[2][1]
    def cmp(a,b):
      if a[1] < b[1]: return 1
      elif a[1] > b[1]: return -1
      return 0
    n = sorted(n, cmp)
    m = None
    original = (n[1][1], n[2][1])
    max0 = 1+max( int( n[1][1] * floor(n[0][1]/n[1][1] + 1) ), 
                  int( n[2][1] * floor(n[0][1]/n[2][1] + 1) ) )
    if n[0][1] == n[1][1]: max0 = 1 + int( n[2][1] * floor(n[0][1]/n[2][1] + 1) )

#   max0 = n[0][1]+1
#   n[2][1] = n[0][1]
#   n[1][1] = n[0][1]
    for n[0][1] in range(n[0][1], max0): # not sure any of these loops help
      for n[1][1] in range(n[1][1], n[0][1]+1):
        for n[2][1] in range(n[2][1], n[1][1]+1):
          # create grid
          grid = cell.copy()
          for i in range(3): grid[:,n[i][0]] *= 1e0 / float(n[i][1] )
          if not Density._goodsyms(grid, pointgroup): continue

          if m ==  None: m = grid.copy()
          if np.linalg.det(grid) < np.linalg.det(m):  m = grid.copy()

        n[2][1] = original[1]
      n[1][1] = original[0]
        
    assert m != None

    grid = m
    for i in range(3): result += "%18.12f %18.12f %18.12f\n" % tuple(grid[:,i])
    result += "  %18.12f %18.12f %18.12f\n" % tuple(self.offset)
    return result


  @staticmethod
  def _goodsyms(grid, pointgroup):
    """ Checks if grid has all symmetries of pointgroup """
    import numpy as np

    def is_unimodular(matrix):
      from math import fabs, floor

      d = np.linalg.det(matrix)
      if not (fabs(d-1e0) < 1e-12 or fabs(d+1e0) < 1e-12): return False
      for i in range(3):
        for i in range(3):
          if fabs( matrix[i,j] - float(floor(matrix[i,j]+0.1)) ) > 1e-12: 
            return False
      return True


    for transform in pointgroup:
      assert np.linalg.norm(transform.trans) < 1e-12
      op = np.matrix( [[transform.op[j,i] for j in range(3)] for i in range(3)] )
      if not is_unimodular( grid.I * op * grid ): return False
    return True

   



