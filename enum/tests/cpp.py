def test_isinteger():
  from numpy import array
  from lada.enum.cppwrappers import is_integer
  assert not is_integer(array([0.5, 0.5]))
  assert not is_integer(array([1.0, 0.5]))
  assert not is_integer(array(0.1))
  assert is_integer(array([1.0, -1.0])) 
  assert is_integer(array(-5))
  assert is_integer(array([ [1.0, -1.0], [5.0, 1111.0] ]))

if __name__ == '__main__':
  test_isinteger()
