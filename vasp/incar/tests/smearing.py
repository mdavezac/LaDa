""" Check smearing special parameter. """
def test():
  from quantities import eV, J
  from pickle import loads, dumps
  from pylada.vasp.incar import Smearing

  assert Smearing(None).value is None
  assert Smearing("fermi").value == "fermi"
  assert Smearing("gaussian").value == "gaussian"
  assert Smearing("mp 1").value == "metal"
  assert Smearing("mp    1").value == "metal"
  assert Smearing("mp1").value == "metal"
  assert Smearing("mp 2").value == "mp 2"
  assert Smearing("mp 3").value == "mp 3"
  assert Smearing("mp   3").value == "mp 3"
  assert Smearing("mp3").value == "mp 3"
  assert Smearing("bloechl").value == "insulator"
  assert Smearing("insulator").value == "insulator"
  assert Smearing("tetra").value == "tetra"
  assert Smearing("metal", 0.2*eV).value[0] == 'metal'\
         and abs(Smearing("metal", 0.2*eV).value[1] - 0.2*eV) < 1e-8
  assert Smearing("metal", 0.3*eV).value[0] == 'metal' \
         and abs(Smearing("metal", 0.3*eV).value[1] - 0.3*eV) < 1e-8
  assert Smearing("metal", (0.3*eV).rescale(J)).value[0] == 'metal' \
         and abs((Smearing("metal", (0.3*eV).rescale(J)).value[1] - 0.3*eV).rescale(eV)) < 1e-8
  assert repr(Smearing("metal", (0.3*eV).rescale(J))) \
           == "Smearing('metal', array(4.80652959e-20) * J)"
  assert Smearing(None, 0.3*eV).value[0] is None \
         and abs(Smearing(None, 0.3*eV).value[1] - 0.3*eV) < 1e-8

  assert Smearing().incar_string() is None
  assert Smearing(None).incar_string() is None
  assert Smearing('metal', (0.3*eV).rescale(J)).incar_string()\
           == 'ISMEAR = 1\nSIGMA = 0.3'
  assert Smearing('fermi').incar_string() == 'ISMEAR = -1'
  assert Smearing('gaussian').incar_string() == 'ISMEAR = 0'
  assert Smearing('bloechl').incar_string() == 'ISMEAR = -5'
  assert Smearing('dynamic').incar_string() == 'ISMEAR = -3'
  assert Smearing('tetra').incar_string() == 'ISMEAR = -4'
  assert Smearing('insulator').incar_string() == 'ISMEAR = -5'
  assert Smearing('mp 1').incar_string() == 'ISMEAR = 1'
  assert Smearing('mp 2').incar_string() == 'ISMEAR = 2'
  assert Smearing('mp 3').incar_string() == 'ISMEAR = 3'
  try: Smearing('mp 4').incar_string()
  except: pass
  else: raise RuntimeError()
  try: Smearing(-6).incar_string()
  except: pass
  else: raise RuntimeError()
  try: Smearing(4).incar_string()
  except: pass
  else: raise RuntimeError()
  assert Smearing(-5).incar_string() == 'ISMEAR = -5'
  assert Smearing(-3).incar_string() == 'ISMEAR = -3'
  assert Smearing(-4).incar_string() == 'ISMEAR = -4'
  assert Smearing(-1).incar_string() == 'ISMEAR = -1'
  assert Smearing(0).incar_string() == 'ISMEAR = 0'
  assert Smearing('mp 1').incar_string() == 'ISMEAR = 1'
  assert Smearing('mp 2').incar_string() == 'ISMEAR = 2'
  assert Smearing('mp 3').incar_string() == 'ISMEAR = 3'

  assert Smearing('dynamic', [0.3, 0.5, 0.6, 0.8]).incar_string()\
           == 'ISMEAR = -3\nSIGMA = 0.3 0.5 0.6 0.8'
  assert Smearing(sigma=[0.3, 0.5, 0.6, 0.8]).incar_string()\
           == 'ISMEAR = -3\nSIGMA = 0.3 0.5 0.6 0.8'
  try: Smearing(0, [0.3, 0.5, 0.6, 0.8])
  except: pass
  else: raise RuntimeError()
  assert loads(dumps(Smearing(sigma=[0.3, 0.5, 0.6, 0.8]))).incar_string()\
           == 'ISMEAR = -3\nSIGMA = 0.3 0.5 0.6 0.8'
  

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

