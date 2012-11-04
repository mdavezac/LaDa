def test_smith():
  from numpy import dot, all
  from numpy.random import randint
  from numpy.linalg import det
  from lada.math import smith_normal_form
  from lada.error import input
  
  for i in xrange(50):
    cell = randint(-5, 5, size=(3,3))
    while det(cell) == 0: cell = randint(-5, 5, size=(3,3))
    s, l, r = smith_normal_form(cell)
    assert all(dot(dot(l, cell), r) == s)

  try: smith_normal_form([[0, 0, 0], [1, 2, 0], [3, 4, 5]])
  except input: pass
  else: raise Exception()

if __name__ == '__main__':
  test_smith()
