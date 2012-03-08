""" Check relaxation property. """
def test():
  from lada.vasp.incar import Relaxation
  from collections import namedtuple
  Vasp = namedtuple('Vasp', ['ediff', 'ediffg'])
  vasp = Vasp(1e-5, 1e-4)

  assert Relaxation(None).value is None
  assert Relaxation("static").value    == "static"
  assert Relaxation("cellshape").value == "cellshape"
  assert Relaxation("ionic").value     == "ionic"
  assert Relaxation("volume").value    == "volume"
  assert set(Relaxation("ionic cellshape").value.split()) == set(["cellshape", 'ionic'])
  assert set(Relaxation("ionic cellshape volume").value.split()) == set(["cellshape", 'ionic', 'volume'])
  try: set(Relaxation("ionic volume").value.split())
  except: pass
  else: raise RuntimeError()

  assert Relaxation('ionic', 50).value == 'ionic'
  assert Relaxation('ionic', 60).value == ('ionic', 60)
  assert Relaxation('ionic cellshape', 60, -1).value == ('static', 60)
  assert Relaxation('ionic cellshape', 60, 2).value == ('ionic cellshape', 60, 2)
  assert Relaxation('cellshape ionic', 60, 2).value == ('ionic cellshape', 60, 2)
  assert Relaxation('cellshape volume ionic', 60, 2, 50).value == ('ionic cellshape volume', 60, 2, 50)
  assert repr(Relaxation('cellshape volume ionic', 60, 2, 50))\
            == "Relaxation('ionic cellshape volume', 60, 2, 50.0)"

  assert Relaxation(None).incar_string(vasp=vasp, structure=range(5)) is None
  assert Relaxation('static').incar_string(vasp=vasp, structure=range(5))\
           == 'ISIF = 2\nIBRION = -1'
  assert Relaxation('ionic', 60).incar_string(vasp=vasp, structure=range(5))\
           == 'ISIF = 2\nNSW = 60'
  assert Relaxation('cellshape', 60).incar_string(vasp=vasp, structure=range(5))\
           == 'ISIF = 5\nNSW = 60'
  assert Relaxation('cellshape ionic', 60).incar_string(vasp=vasp, structure=range(5))\
           == 'ISIF = 4\nNSW = 60'
  assert Relaxation('cellshape ionic volume', 60, 3).incar_string(vasp=vasp, structure=range(5))\
           == 'ISIF = 3\nNSW = 60\nIBRION = 3'
  assert Relaxation('cellshape ionic volume', 60, 3, 20).incar_string(vasp=vasp, structure=range(5))\
           == 'ISIF = 3\nNSW = 60\nPOTIM = 20.0\nIBRION = 3'

  vasp = Vasp(1e-4, 1e-5)
  assert Relaxation(None).incar_string(vasp=vasp, structure=range(5)) is None
  assert Relaxation('static').incar_string(vasp=vasp, structure=range(5))\
           == 'ISIF = 2\nIBRION = -1'
  try: Relaxation('ionic', 60).incar_string(vasp=vasp, structure=range(5)) 
  except: pass
  else: raise RuntimeError()
  vasp = Vasp(1e-4, -1e-5)
  assert Relaxation('ionic', 60).incar_string(vasp=vasp, structure=range(5)) == 'ISIF = 2\nNSW = 60'
  vasp = Vasp(-1e-4, 1e-5)
  try: Relaxation('ionic', 60).incar_string(vasp=vasp, structure=range(5))
  except: pass
  else: raise RuntimeError()
  vasp = Vasp(-1e-6, 1e-4)
  assert Relaxation('cellshape ionic volume', 60, 3, 20).incar_string(vasp=vasp, structure=range(5))\
           == 'ISIF = 3\nNSW = 60\nPOTIM = 20.0\nIBRION = 3'


if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

