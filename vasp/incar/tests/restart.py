def test():
  from collections import namedtuple
  from pickle import loads, dumps
  from lada.vasp.incar._params import Restart

  Extract = namedtuple("Extract", ['directory'])
  Vasp = namedtuple("Vasp", ['nonscf'])

  assert Restart(None).incar_string(vasp=None) \
            ==  "ISTART = 0   # start from scratch\n"\
                "ICHARG = 2   # superpositions of atomic densities"
  assert Restart(None).incar_string(vasp=Vasp(True)) \
            ==  "ISTART = 0   # start from scratch\n"\
                "ICHARG = 12   # superpositions of atomic densities"
  assert Restart(None).incar_string(vasp=Vasp(False)) \
            ==  "ISTART = 0   # start from scratch\n"\
                "ICHARG = 2   # superpositions of atomic densities"





if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

