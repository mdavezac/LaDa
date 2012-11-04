def test_gruber():
  from numpy import dot
  from numpy.linalg import inv
  from lada.math import gruber, is_integer
  from lada.error import internal, input
  
  cell = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
  lim = 5
   
  for a00 in [-1, 1]: 
    for a10 in xrange(-lim, lim+1): 
      for a11 in [-1, 1]:
        for a20 in xrange(-lim, lim+1): 
          for a21 in xrange(-lim, lim+1): 
            for a22 in [-1, 1]:
              a = [[a00, 0, 0], [a10, a11, 0], [a20, a21, a22]]
              g = gruber(dot(cell, a))
              assert is_integer(dot(inv(cell), g))
              assert is_integer(dot(inv(g), cell))

  try: gruber([[0, 0, 0], [1, 2, 0], [4, 5, 6]])
  except input: pass
  else: raise Exception()

  try: gruber([[1, 0, 0], [1, 1, 0], [4, 5, 1]], itermax=2)
  except internal: pass
  else: raise Exception()

if __name__ == '__main__':
  test_gruber()
