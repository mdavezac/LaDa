def test():
  from collections import namedtuple
  from pickle import loads, dumps
  from lada.vasp.incar._params import Npar

  Comm = namedtuple('Comm', ['n'])

  # ispin == 1
  assert Npar(0).incar_string(comm=Comm(16)) is None
  assert Npar(1).incar_string(comm=Comm(16)) == "NPAR = 1"
  assert Npar(2).incar_string(comm=Comm(16)) == "NPAR = 2"
  assert Npar('power of two').incar_string(comm=Comm(16)) == "NPAR = 4"
  assert Npar('power of two').incar_string(comm=Comm(17)) is None
  assert Npar('power of two').incar_string(comm=Comm(2)) == "NPAR = 1"
  assert Npar('sqrt').incar_string(comm=Comm(16)) == "NPAR = 4"
  assert repr(Npar(1)) == "Npar(1)"
  assert repr(loads(dumps(Npar(2)))) == "Npar(2)"
  assert repr(loads(dumps(Npar("power of two")))) == "Npar('power of two')"

if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

