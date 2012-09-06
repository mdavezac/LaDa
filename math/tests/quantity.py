""" Checks atom methods and attributes. """
def test():
  """ Test atom initialization. """
  from _quantity import get_static_object, check_get, giveinRy
  from numpy import abs
  from quantities import meter, cm, eV

  a = get_static_object()
  assert a.__class__ is (1*meter).__class__
  assert abs(a - (1*meter)) < 1e-12

  assert abs(giveinRy(1.0) - 1e0) < 1e-8
  assert abs(giveinRy(1L) - 1e0) < 1e-8
  assert abs(giveinRy(2*eV) - (2*eV).rescale("Ry").magnitude.tolist()) < 1e-8


if __name__ == "__main__":
  from lada.crystal.cppwrappers import Atom
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  test()
