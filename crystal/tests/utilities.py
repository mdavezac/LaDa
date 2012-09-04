""" Test some utilities from crystal extension module. """
def test_periodic():
  """ Test periodic images. """
  from lada.crystal.cppwrappers import are_periodic_images
  from numpy import array, dot
  from numpy.linalg import inv
  from random import randint, random
  cell = array([ [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0] ])
  invcell = inv(cell)
  for i in range(10):
    vec = dot(cell, array([random(), random(), random()]))
    for j in range(10):
      vec2 = vec + dot(cell, array([randint(-10, 11), randint(-10, 11), randint(-10, 11)]))
      assert are_periodic_images(vec, vec2, invcell)
      vec3 = vec2 + dot(cell, array([random()+0.0001, random(), random()]))
      assert not are_periodic_images(vec, vec3, invcell)

def test_into_cell():
  """ Test that vector is folded to fractional coordinates >= 0 and < 1. """
  from lada.crystal.cppwrappers import are_periodic_images, into_cell
  from numpy import array, dot, all
  from numpy.linalg import inv
  from random import uniform
  cell = array([ [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0] ])
  invcell = inv(cell)
  for i in range(100):
    vec = dot(cell, array([uniform(-10, 10), uniform(-10, 10), uniform(-10, 10)]))
    vec2 = into_cell(vec, cell, invcell)
    assert are_periodic_images(vec, vec2, invcell)
    assert all(dot(invcell, vec2) >= 0e0) and all(dot(invcell, vec2) < 1e0)
  
def test_zero_centered():
  """ Test that vector is folded to fractional coordinates >= -0.5 and < 0.5. """
  from lada.crystal.cppwrappers import are_periodic_images, zero_centered
  from numpy import array, dot, all
  from numpy.linalg import inv
  from random import uniform
  cell = array([ [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0] ])
  invcell = inv(cell)
  for i in range(100):
    vec = dot(cell, array([uniform(-10, 10), uniform(-10, 10), uniform(-10, 10)]))
    vec2 = zero_centered(vec, cell, invcell)
    assert are_periodic_images(vec, vec2, invcell)
    assert all(dot(invcell, vec2) >= -0.5) and all(dot(invcell, vec2) < 0.5)

def test_into_voronoi():
  """ Test that vector is folded into first Brillouin zone. """
  from lada.crystal.cppwrappers import are_periodic_images, into_voronoi
  from numpy import array, dot
  from numpy.linalg import inv
  from random import uniform
  cell = array([ [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0] ])
  invcell = inv(cell)
  for i in range(100):
    vec = dot(cell, array([uniform(-10, 10), uniform(-10, 10), uniform(-10, 10)]))
    vec2 = into_voronoi(vec, cell, invcell)
    assert are_periodic_images(vec, vec2, invcell)
    n = dot(vec2, vec2)
    assert dot(vec, vec) >= n
    for j in range(-3, 4):
      for k in range(-3, 4):
        for l in range(-3, 4):
          o = vec2 + dot(cell, [j, k, l])
          assert dot(o, o) >= n

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])

  test_periodic()
  test_into_cell()
  test_zero_centered()
  test_into_voronoi()
