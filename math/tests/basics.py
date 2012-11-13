def test_isinteger():
  from numpy import array
  from lada.math import is_integer
  from lada.error import TypeError

  assert is_integer(array([0, 1, 0, 2]))
  assert is_integer(array([0L, 1, 0, 2]))
  assert is_integer(array([0, 1.0, 0.0, 2.0]))
  assert not is_integer(array([0, 1.1, 0.0, 2.0]))
  try: is_integer([5, 6, 7])
  except TypeError: pass
  else: raise Exception()

def test_floorint():
  from numpy import array, all 
  from lada.math import floor_int
  from lada.error import TypeError

  assert floor_int(array([0.1, -0.1, -0.5, 0.5, -0.55, 0.55])).dtype == 'int64'
  assert all( floor_int(array([0.1, -0.1, -0.5, 0.5, -0.55, 0.55, 1.999, -1.99]))
               == [0, -1, -1, 0, -1, 0, 1, -2] )
  assert all( floor_int(array([[0.1, -0.1, -0.5, 0.5], [-0.55, 0.55, 1.999, -1.99]]))
               == [[0, -1, -1, 0], [-1, 0, 1, -2]] )
  
  try: floor_int([5, 6, 7])
  except TypeError: pass
  else: raise Exception()

def test_rotation():
  from numpy import all, abs, pi
  from lada.math import Rotation
  from lada.error import TypeError

  assert all(abs(Rotation(0, [1, 0, 0]) - [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]]) < 1e-8)
  assert all(abs(Rotation(pi, [1, 0, 0]) - [[1, 0, 0], [0, -1, 0], [0, 0, -1], [0, 0, 0]]) < 1e-8)
  assert all(abs(Rotation(pi/2, [1, 0, 0]) - [[1, 0, 0], [0, 0, -1], [0, 1, 0], [0, 0, 0]]) < 1e-8)

  try: Rotation(pi/2, [0, 0])
  except TypeError: pass
  else: raise Exception()
  try: Rotation(pi/2, [0, 0, 1, 0])
  except TypeError: pass
  else: raise Exception()

def test_translation():
  from numpy import all, abs, identity, pi
  from lada.math import Translation
  from lada.error import TypeError

  assert all(abs(Translation([2, 2, 2])[:3] - identity(3)) < 1e-8)
  assert all(abs(Translation([2, 2, 2])[3] - 2) < 1e-8)
  assert all(abs(Translation([pi, pi/2., 2])[:3] - identity(3)) < 1e-8)
  assert all(abs(Translation([pi, pi/2., 2])[3] - [pi, pi/2., 2]) < 1e-8)

  try: Translation([0, 0])
  except TypeError: pass
  else: raise Exception()
  try: Translation([0, 0, 1, 0])
  except TypeError: pass
  else: raise Exception()

if __name__ == '__main__': 
  test_isinteger()
  test_floorint()
  test_rotation()
  test_translation()
