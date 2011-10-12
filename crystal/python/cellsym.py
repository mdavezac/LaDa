#! /usr/bin/python

def lattice_type( cell, tiny = 1e-8 ):
  """ Returns lattice type  for specifically parameterized cells 
      The return consists of the 3-tuple (bravais lattice, celldimensions, new cell),
      where bravais lattice is an integer between 1 and 15, 
      the cell dimensions are (a, b/a, c/a, cos(alpha), cos(beta), cos(gamma)),
      and the new cell may have had its vectors swapped. 
  """
  from numpy import matrix, zeros
  from numpy.linalg import norm
  from math import fabs, sqrt

  matcell = matrix( cell.copy() )
  lengths = matcell.T * matcell 

  lengths[0,1] = lengths[0,1] / lengths[0,0] / lengths[1,1]
  lengths[1,0] = lengths[0,1]
  lengths[0,2] = lengths[0,2] / lengths[0,0] / lengths[2,2]
  lengths[2,0] = lengths[2,0]
  lengths[1,2] = lengths[1,2] / lengths[1,1] / lengths[2,2]
  lengths[2,1] = lengths[1,2]

  a1 = matcell[:,0] # sometimes more expressive to recompute
  a2 = matcell[:,1]
  a3 = matcell[:,2]

  def is_bodycentered():
    return fabs(norm(a1) - norm(a2)) < tiny \
       and fabs(norm(a2) - norm(a3)) < tiny \
       and fabs((a1+a2) * (a1*a3)) < tiny \
       and fabs((a1+a2) * (a2*a3)) < tiny \
       and fabs((a1+a3) * (a2*a3)) < tiny \
       and fabs(norm(a1+a3) - norm(a2*a3)) < tiny \
       and norm(a1*a2) > norm(a1+a3) + tiny

  def is_facecentered():
    return fabs(norm(a1) - norm(a2-a3)) < tiny \
       and fabs(norm(a2) - norm(a1-a3)) < tiny \
       and fabs(norm(a3) - norm(a1-a2)) < tiny \
       and norm(a1+a2-a3) > norm(a1+a3-a2) + tiny \
       and norm(a1+a3-a2) > norm(a2+a3-a1) + tiny  

  def cosine(a, b): return a.T * b / norm(a) / norm(b)

  ibrav = 15
  celldims = zeros( (6,), dtype="float64" )
  print cosine(a1, a2), cosine(a1, a3), cosine(a2, a3), norm(a1), norm(a2)
  if fabs(cosine(a1, a2) - cosine(a1, a3)) < tiny  and fabs(cosine(a1, a2) - cosine(a2, a3)) < tiny:
    # a1 = a2 = a3
    if fabs(norm(a1)-norm(a2)) < tiny and fabs(norm(a2)-norm(a3)) < tiny:
      if fabs(cosine(a1, a2)) < tiny:
        ibrav = 1
        celldims[0] = norm(a1)
      elif fabs(cosine(a1, a2) + 0.33e0) < tiny: 
        ibrav = 2
        celldims[0] = norm(a1) * 2e0 / sqrt(3e0)
      elif fabs(cosine(a1, a2) - 0.5e0) < tiny: 
        ibrav = 3
        celldims[0] = norm(a1) * sqrt(2e0)
      else: 
        ibrav = 7
        celldims[0] = norm(a1)
        celldims[3] = cosine(a1, a2)
    # a1 perpendicular to a2 (hence a2 to a3 according to first condition)
    elif fabs(cosine(a1, a2)) < tiny:
      if fabs(norm(a1)-norm(a2)) < tiny:
        ibrav = 5
        celldims[0] = norm(a1)
        celldims[2] = norm(a3) / norm(a1)
      elif norm(a3) > norm(a2) + tiny and norm(a2) > norm(a1) + tiny:
        ibrav = 8
        celldims[0] = norm(a1)
        celldims[1] = norm(a2) / norm(a1)
        celldims[2] = norm(a3) / norm(a1) 
  elif fabs(cosine(a3, a1)-cosine(a3, a2)) < tiny:
    # one othogonal angle
    if fabs(cosine(a3, a1)) < tiny: 
      # a = b
      if fabs(norm(a1) - norm(a2)) < tiny:
        # Cosine(alpha) equal to -1/2 means hexagonal: (IBRAV=4)
        # Adjustment: 'c-axis' shall be the special axis.
        if fabs(cosine(a2, a1) + 0.5e0) < tiny:
           ibrav=4
           celldims[0] = norm(a1)
           celldims[2] = norm(a3) / norm(a1) 
        # Other angles mean base-centered orthorhombic: (IBRAV=11)
        # Adjustment: Cosine between A1 and A2 shall be lower than zero, the
        #             'c-axis' shall be the special axis.
        elif cosine(a2, a1) < -tiny:
          ibrav=11
          celldims[0] = norm(a2+a1)
          celldims[1] = norm(a2-a1) / celldims[0]
          celldims[2] = norm(a3) / celldims[0]
      elif cosine(a2, a1) < -tiny and norm(a1) > norm(a2) + tiny:
        ibrav = 12
        print "Warning: Monoclinic adjustments (a1->a3, a2->a1, a3->a2).\n"
        celldims[0] = norm(a2)
        celldims[1] = norm(a3) / celldims[0]
        celldims[2] = norm(a1) / celldims[0]
        celldims[5] = lengths[0,1]
        old_cell = matcell.copy()
        matcell[:,2] = old_cell[:,0]
        matcell[:,0] = old_cell[:,1]
        matcell[:,1] = old_cell[:,3]
    else: # arbitrary angles
      if is_bodycentered():
        ibrav = 6
        celldims[0] = norm(a3+a1)
        celldims[2] = norm(a2+a1) / celldims[0]
      elif fabs(norm(a1) - norm(a2)) < tiny and cosine(a1, a3) < -tiny and cosine(a2, a3) < -tiny:
        ibrav = 13
        celldims[0] = norm(a2+a1)
        celldims[1] = norm(a2-a1) / celldims[0]
        celldims[2] = norm(a3) / celldims[0]
        celldims[5] = (a3.T*a2 + a3.T*a1) / norm(a2+a1) / norm(a3) 
  else:
    if is_bodycentered():
      ibrav = 9
      celldims[0] = norm(a3+a2)
      celldims[2] = norm(a3+a1) / celldims[0]
      celldims[3] = norm(a2+a1) / celldims[0]
    if is_facecentered():
      ibrav = 10
      celldim[0] = norm(a2+a3-a1)
      celldim[1] = norm(a3+a1-a2) / celldim[0]
      celldim[2] = norm(a1+a2-a3) / celldim[0]
      
  return ibrav, celldims, matcell



def lattice_symmetry(cell, tiny = 1e-8):
  """ Computes the cell with smallest possible column-vectors, and returns lattice-type. 
      First looks for smallest possible vectors.
      Then performs all type of rotations to find the cell parameterized as expected by 
      latice_type() subroutine. 
      Return is that of lattice_type(), where the cell may be a linear
      combination and rotation of the input cell-vectors.
  """

  from math import fabs
  from numpy.linalg import norm, det

  # make cell direct:
  def direct(cell):
    from numpy.linalg import det 
    from math import pow
    if det(cell) < 0e0: cell = -cell
    assert pow( det(cell), 0.33333 ) > tiny, "Cell has zero volume."
    return cell

  cell = direct(cell)

  # get lattice type and dimensions. 
  ibrav_in, celldims_in, cell_in = lattice_type(cell, tiny)

  a1 = cell_in[:,0]
  a2 = cell_in[:,1]
  a3 = cell_in[:,2]
  print "in: ", ibrav_in
  print cell_in
  print det(cell_in)


  # looks for equivalent cell with shortest cell-vectors.
  def try_minus(a, b, changed = False, first_try = True):
    """ Tries to increase or decrease a by b until shortest lenght is found.
        If fails on first try, then tries to increase a by b.
        returns new a, and whether there has been any change.
    """
    n1, n0 = norm(a-b), norm(a)
    if n1 >= n0  + tiny:
      n1, n0 = norm(a+b), norm(a)
      while n1 < n0 + tiny:
        a += b
        if n1 < n0 - tiny: changed = True
        n1, n0 = norm(a+b), norm(a)
      return a, changed
    
    while n1 < n0 + tiny:
      a -= b
      if n1 < n0 - tiny: changed = True
      n1, n0 = norm(a-b), norm(a)
    return a, changed

  changed = True
  while( changed ):
    a1, changed = try_minus(a1, a2)
    a1, changed = try_minus(a1, a3, changed)
    a2, changed = try_minus(a2, a1, changed)
    a2, changed = try_minus(a2, a3, changed)
    a3, changed = try_minus(a3, a1, changed)
    a3, changed = try_minus(a3, a2, changed)

  # now tries transformations.
  def unimodular_transform(transformation = None, axis = None):
    """ Generator to loop through cell with integer elements in [-2,2].
        Avoids matrices whith det() != 1
    """
    from numpy.linalg import det
    from numpy import zeros, matrix

    if transformation is None: transformation = matrix(zeros((3,3),dtype="float64"))
    if axis is None: axis = (2,2) 
    
    next_axis = None
    if axis[1] != -2: next_axis = (axis[0], axis[1]-1)
    elif axis[0] != -2: next_axis = (axis[0]-1, 0)

    for i in range(-2,3):
      transformation[axis] = i
      if next_axis is not None:
        for transformation in unimodular_transform(transformation, next_axis):
          if det(transformation) == 1: yield transformation
      elif det(transformation) == 1: yield transformation

  ibrav, celldims, new_cell = lattice_type(cell, tiny)
  for transform in unimodular_transform():
    ibrav_try, celldims_try, new_cell_try = lattice_type(cell*transform, tiny)
    if    ibrav_try < ibrav \
       or (     ibrav_try == ibrav
            and fabs(celldims_try[3]) < fabs(celldims[3]) - tiny 
            and fabs(celldims_try[4]) < fabs(celldims[4]) - tiny 
            and fabs(celldims_try[5]) < fabs(celldims[5]) - tiny ):
      print ibrav_try
      print transform
      ibrav = ibrav_try
      celldims = celldims_try.copy()
      new_cell = new_cell_try.copy()

  # Finally keeps input answer if equivalent to this answer.
  if ibrav_in == ibrav:
    celldims = celldims_in.copy()
    new_cell = cell_in.copy()

  return ibrav, celldims, new_cell


def str_lattice_type(ibrav, celldims):
  """ Prints the lattice type in human readable format. """

  if   ibrav ==  1:
    return "Cubic Cell (a=%f)" % (celldims[0])
  elif ibrav ==  2:
    return "Body-centered cubic cell (a=%f)" % (celldims[0])
  elif ibrav ==  3:
    return "Face-centered cubic cell (a=%f)" % (celldims[0])
  elif ibrav ==  4:
    return "Hexagonal cell (a=%f, a/c=%f)" % (celldims[0], celldims[2])
  elif ibrav ==  5:
    return "Simple tetragonal cell (a=%f, a/c=%f)" % (celldims[0], celldims[2])
  elif ibrav ==  6:
    return "Body-centered tetragonal cell (a=%f, a/c=%f)" % (celldims[0], celldims[2])
  elif ibrav ==  7:
    return "Romboedric cell (a=%f, cos(alpha)=%f)" % (celldims[0], celldims[3])
  elif ibrav ==  8:
    return "Simple orthorombic cell (a=%f, b/a=%f, c/a=%f)" \
           % (celldims[0], celldims[1], celldims[2])
  elif ibrav ==  9:
    return "Body-centered orthorombic cell (a=%f, b/a=%f, c/a=%f)" \
           % (celldims[0], celldims[1], celldims[2])
  elif ibrav == 10:
    return "Face-centered orthorombic cell (a=%f, b/a=%f, c/a=%f)" \
           % (celldims[0], celldims[1], celldims[2])
  elif ibrav == 11:
    return "Base-centered orthorombic cell (a=%f, b/a=%f, c/a=%f)" \
           % (celldims[0], celldims[1], celldims[2])
  elif ibrav == 12:
    return "Simple monoclinic cell (a=%f, b/a=%f, c/a=%f, cos(beta)=%f)" \
           % (celldims[0], celldims[1], celldims[2], celldims[4])
  elif ibrav == 13:
    return "Base-centered monoclinic cell (a=%f, b/a=%f, c/a=%f, cos(beta)=%f)" \
           % (celldims[0], celldims[1], celldims[2], celldims[4])
  elif ibrav == 14:
    return "Triclinic cell (a=%f, b/a=%f, c/a=%f, cos(alpha)=%f, cos(beta)=%f, cos(gamma)=%f)" \
           % (celldims[0], celldims[1], celldims[2], celldims[3], celldims[4], celldims[5])
  elif ibrav == 15: return "Unkown ibrav type: %i.\n"  % (ibrav)

    

def main():
  from numpy import array
  from numpy import linalg


  cell = matrix(array( [ [0.0, 8, 8], [1, 7, 8], [1, 1, 0] ] ) * 4.419999999999999929)
  w, v = linalg.eig(cell)
  print w, v

# cell = array( [ [0.0, 1, 1], [1, 0, 1], [1, 1, 0] ] )
# structure = crystal.read_poscar( ("K", "Rb"), "POSCAR" )
# print structure.cell, "\n"
# ibrav, celldims, cell = lattice_symmetry( cell )
# print str_lattice_type( ibrav, celldims )



if __name__ == "__main__" :
  main()
