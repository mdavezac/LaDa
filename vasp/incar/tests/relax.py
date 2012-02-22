""" Check relaxation property. """
def test():
  from lada.vasp.incar import Incar

  a = Incar()
  assert a.relaxation == "static"
  assert a.isif == 1 and (a.ibrion == -1 or a.ibrion is None) and (a.nsw is None or a.nsw == 0)
  assert a.relaxation == "static"
  a.relaxation = "volume"
  assert a.relaxation == "volume"
  assert a.isif == 7 and a.ibrion == 2 and a.nsw == 50
  a.relaxation = "ions"
  assert a.relaxation == "ionic"
  assert a.isif == 1 and a.ibrion == 2 and a.nsw == 50
  a.relaxation = "cellshape"
  assert a.relaxation == "cellshape"
  assert a.isif == 5 and a.ibrion == 2 and a.nsw == 50
  a.relaxation = "volume"
  assert a.relaxation == "volume"
  assert a.isif == 7 and a.ibrion == 2 and a.nsw == 50
  a.relaxation = "volume cellshape"
  assert set(a.relaxation.split()) == set(["volume", "cellshape"])
  assert a.isif == 6 and a.ibrion == 2 and a.nsw == 50
  a.relaxation = "ionic cellshape"
  assert set(a.relaxation.split()) == set(["ionic", "cellshape"])
  assert a.isif == 4 and a.ibrion == 2 and a.nsw == 50
  a.relaxation = "ionic cellshape volume"
  assert set(a.relaxation.split()) == set(["ionic", "cellshape", "volume"])
  assert a.isif == 3 and a.ibrion == 2 and a.nsw == 50
  try: a.relaxation = "ionic volume"
  except: pass
  else: raise RuntimeError()
  a.relaxation = "ionic cellshape volume", 3
  try: a.relaxation = "wtf"
  except: pass
  else: raise RuntimeError()

if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

