""" Check relaxation property. """
def test():
  from lada.vasp.incar import Relaxation

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

  assert Relaxation(('ionic', 50)).value == 'ionic'
  assert Relaxation(('ionic', 60)).value == ('ionic', 60)
  assert Relaxation(('ionic cellshape', 60, -1)).value == ('static', 60)
  assert Relaxation(('ionic cellshape', 60, 2)).value == ('ionic cellshape', 60, 2)
  assert Relaxation(('cellshape ionic', 60, 2)).value == ('ionic cellshape', 60, 2)
  assert Relaxation(('cellshape volume ionic', 60, 2, 50)).value == ('ionic cellshape volume', 60, 2, 50)
  assert repr(Relaxation(('cellshape volume ionic', 60, 2, 50)))\
            == "Relaxation(('ionic cellshape volume', 60, 2, 50))"


if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

