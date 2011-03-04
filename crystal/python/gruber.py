""" Contains a class to determine the Niggli cell of a structure.
    Adapted from cctbx/uctbx/gruber_1973.py
"""


class Reduction(object):
  """ Returns Niggli unit cell. """

  def __init__(self, itermax=None, multiplier=16, nochange=2):
    from numpy import matrix

    self.itermax = itermax
    if (itermax is None): self.itermax = 100

    self.multiplier = multiplier
    self.nochange   = nochange

  def __call__(self, cell, recip=False):
    """ Returns Niggli unit cell. """

    from copy import deepcopy
    from numpy import matrix
    from numpy.linalg import inv

    cell = matrix(cell).copy()
    if recip: cell = inv(cell).T

    metrical = cell.T * cell

    self.a = metrical[0,0]
    self.b = metrical[1,1]
    self.c = metrical[2,2]
    self.d = 2e0 * metrical[1,2]
    self.e = 2e0 * metrical[0,2]
    self.f = 2e0 * metrical[0,1]
    self._r_inv = matrix( [[1,0,0],[0,1,0],[0,0,1]], dtype="float64")
    self._n_iterations = 0
    self._n_no_significant_change = 0
    self._last_abc_significant_change_test = (-self.a,-self.b,-self.c)

    while (self._step()): pass

    cell = cell * self._r_inv
    if recip: cell = inv(cell).T

    return cell

  def itermax(self): return self.itermax

  def n_iterations(self): return self._n_iterations

  def _def_test(self):
    n_zero = 0
    n_positive = 0
    if (0 < self.d): n_positive += 1
    elif (not self.d < 0): n_zero += 1
    if (0 < self.e): n_positive += 1
    elif (not self.e < 0): n_zero += 1
    if (0 < self.f): n_positive += 1
    elif (not self.f < 0): n_zero += 1
    return n_zero, n_positive

  def _def_gt_0(self):
    n_zero, n_positive = self._def_test()
    return n_positive == 3 or (n_zero == 0 and n_positive == 1)

  def _step(self):
    # N1
    if self.b < self.a: self._n1_action()
    # N2
    if self.c < self.b:
      self._n2_action()
      return True
    # N3
    if self._n3_action(condition=self._def_gt_0()): return False
    if self._b2_action(): return True
    if self._b3_action(): return True
    if self._b4_action(): return True
    if self._b5_action(): return True
    return False

  def _is_changing(self):
    abc = (self.a,self.b,self.c)
    m = self.multiplier
    change = tuple([(new*m+(new-last))-new*m
      for new,last in zip(abc, self._last_abc_significant_change_test)])
    if (change == (0,0,0)):
      self._n_no_significant_change += 1
      if (self._n_no_significant_change >= self.nochange): return False
    else: self._n_no_significant_change = 0
    self._last_abc_significant_change_test = abc
    return True

  def _cb_update(self, m_elems):
    from numpy import matrix

    if (self._n_iterations == self.itermax):
      raise "Gruber iteration limit exceeded (limit=%d)." % self.itermax 
    self._r_inv *= matrix(m_elems)
    self._n_iterations += 1

  def _n1_action(self):
    from numpy import matrix 

    self._cb_update( matrix([[0,-1,0],[-1,0,0],[0,0,-1]], dtype="float64") )
    self.a, self.b = self.b, self.a
    self.d, self.e = self.e, self.d

  def _n2_action(self):
    from numpy import matrix 

    self._cb_update( matrix([[-1,0,0],[0,0,-1],[0,-1,0]], dtype="float64") )
    self.b, self.c = self.c, self.b
    self.e, self.f = self.f, self.e

  def _n3_action(self, condition):
    from math import fabs as abs
    from numpy import diag

    if condition: 
      i,j,k = 1,1,1
      if (self.d < 0): i = -1
      if (self.e < 0): j = -1
      if (self.f < 0): k = -1
      self._cb_update( diag([i,j,k]) )
      self.d = abs(self.d)
      self.e = abs(self.e)
      self.f = abs(self.f)
      return False
    else:
      f = [1,1,1]
      z = -1
      if (0 < self.d): f[0] = -1
      elif (not self.d < 0): z = 0
      if (0 < self.e): f[1] = -1
      elif (not self.e < 0): z = 1
      if (0 < self.f): f[2] = -1
      elif (not self.f < 0): z = 2
      if (f[0]*f[1]*f[2] < 0):
        assert z != -1
        f[z] = -1
      self._cb_update( diag(f) )
      self.d = -abs(self.d)
      self.e = -abs(self.e)
      self.f = -abs(self.f)
      return not self._is_changing()

  def _b2_action(self):
    from math import fabs as abs
    from numpy import matrix

    if (not self.b < abs(self.d)): return False
    j = self._entier((self.d+self.b)/(2*self.b))
    if (j == 0): return False
    self._cb_update( matrix([[1,0,0],[0,1,-j],[0,0,1]], dtype="float64") )
    self.c += j*j*self.b - j*self.d
    self.d -= 2*j*self.b
    self.e -= j*self.f
    assert 0 < self.c
    return True

  def _b3_action(self):
    from math import fabs as abs
    from numpy import matrix

    if (not self.a < abs(self.e)): return False
    j = self._entier((self.e+self.a)/(2*self.a))
    if (j == 0): return False
    self._cb_update( matrix([[1,0,-j],[0,1,0],[0,0,1]], dtype="float64") )
    self.c += j*j*self.a - j*self.e
    self.d -= j*self.f
    self.e -= 2*j*self.a
    assert 0 < self.c
    return True

  def _b4_action(self):
    from math import fabs as abs
    from numpy import matrix

    if (not self.a < abs(self.f)): return False
    j = self._entier((self.f+self.a)/(2*self.a))
    if (j == 0): return False
    self._cb_update( matrix([[1,-j,0],[0,1,0],[0,0,1]], dtype="float64") )
    self.b += j*j*self.a - j*self.f
    self.d -= j*self.e
    self.f -= 2*j*self.a
    assert 0 < self.b
    return True

  def _b5_action(self):
    from math import fabs as abs
    from numpy import matrix

    de = self.d + self.e
    fab = self.f + self.a + self.b
    if (not de+fab < 0): return False
    j = self._entier((de+fab)/(2*fab))
    if (j == 0): return False
    self._cb_update( matrix( [[1,0,-j],[0,1,-j],[0,0,1]], dtype="float64") )
    self.c += j*j*fab-j*de
    self.d -= j*(2*self.b+self.f)
    self.e -= j*(2*self.a+self.f)
    assert 0 < self.c
    return True

  @staticmethod 
  def _entier(x):
    "greatest integer which is not greater than x"
    result = int(x)
    if x-result < 0: result -= 1
    if not (x-result < 1): result += 1 # work around rounding errors
    return result

def main():

  from numpy import matrix
  from numpy.linalg import inv, det 
  import cellsym
  from lada import crystal, eigen

 #file = open("POSCAR", "w")
 #print >>file
 #print >>file, "4.419999999999999929"
 #print >>file, "     0.000000000000000000       1.000000000000000000       1.000000000000000000"
 #print >>file, "     8.000000000000000000       7.000000000000000000       1.000000000000000000"
 #print >>file, "     8.000000000000000000       8.000000000000000000       0.000000000000000000"
 #print >>file, "5 3"
 #print >>file, "Direct"
 #print >>file, "       0.000000000000000000      0.000000000000000000      0.000000000000000000"
 #print >>file, "       0.000000000000000000      0.000000000000000000      0.125000000000000000"
 #print >>file, "       0.000000000000000000      0.000000000000000000      0.375000000000000000"
 #print >>file, "       0.000000000000000000      0.000000000000000000      0.500000000000000000"
 #print >>file, "       0.000000000000000000      0.000000000000000000      0.750000000000000000"
 #print >>file, "       0.000000000000000000      0.000000000000000000      0.250000000000000000"
 #print >>file, "       0.000000000000000000      0.000000000000000000      0.625000000000000000"
 #print >>file, "       0.000000000000000000      0.000000000000000000      0.875000000000000000"
 #file.close()
 #structure = crystal.read_poscar( ("K", "Rb"), "POSCAR" )
 #print structure.cell, "\n"

  cell = matrix( [ [0.0, 8, 8], [1, 7, 8], [1, 1, 0] ], dtype="float32") # * 4.419999999999999929

  reduction = Reduction()

# cell = inv(cell).T
  new_cell = reduction(cell, recip=True)
  print cell, "\n\n", new_cell
  
  ibrav, celldims, dummy = cellsym.lattice_type(cell) 
  print cellsym.str_lattice_type(ibrav, celldims)
  ibrav, celldims, dummy = cellsym.lattice_type(new_cell) 
  print cellsym.str_lattice_type(ibrav, celldims)
 #for i in range(3):
 #  for j in range(3):
 #    structure.cell[i,j] = new_cell[i,j]

 #crystal.print_poscar(structure, ("K", "Rb"))

if __name__ == "__main__":
  main()
