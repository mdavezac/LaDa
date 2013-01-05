""" Checks that point group of cell is determined correctly. """
def test(cell, numops):
  """ Test structure initialization. """
  from numpy import all, abs, dot, array
  from numpy.linalg import inv, det
  from pylada.crystal.cppwrappers import cell_invariants, Structure
  from pylada.math import is_integer

  ops = cell_invariants(cell) 
  if isinstance(cell, Structure): cell = cell.cell.copy()
  assert len(ops) == numops
  for op in ops:
    assert op.shape == (4, 3)
    assert all(abs(op[3, :]) < 1e-8) 
    transformation = dot(dot(inv(cell), op[:3]), cell)
    assert is_integer(transformation)
    assert abs(abs(det(transformation))-1e0) < 1e-8
  if numops != 48:
    allops = cell_invariants(array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]))
    failed = 0
    for op in allops: 
      transformation = dot(dot(inv(cell), op[:3]), cell)
      if not (is_integer(transformation) and abs(abs(det(transformation))-1e0) < 1e-8):
        failed += 1
    assert failed == 48 - numops

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  def test_(cell, numops):
    from pylada.crystal.cppwrappers import Structure
    from numpy import array
    test(array(cell), numops)
    # also good test to make sure that structure storage is correct.
    test(Structure(cell).add_atom(0,0,0, "Si"), numops)

  test_([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]], 48)
  test_([[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]], 48)
  test_([[-0.6,0.5,0.5],[0.6,-0.5,0.5],[0.6,0.5,-0.5]], 4)
  test_([[-0.7,0.7,0.7],[0.6,-0.5,0.5],[0.6,0.5,-0.5]], 8)
  test_([[-0.765,0.7,0.7],[0.665,-0.5,0.5],[0.6,0.5,-0.5]], 2)
