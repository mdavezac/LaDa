""" Check smearing special parameter. """
def test():
  from quantities import eV, J
  from lada.vasp.incar import Smearing

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
  assert Smearing("metal", 0.2*eV).value == 'metal'
  assert Smearing("metal", (0.2*eV).rescale(J)).value == 'metal'
  assert Smearing("metal", 0.3*eV).value[0] == 'metal' \
         and abs(Smearing("metal", 0.3*eV).value[1] - 0.3*eV) < 1e-8
# assert Smearing("metal", (0.3*eV).rescale(J)).value[0] == 'metal' \
#        and abs(Smearing("metal", 0.3*eV).value[1] - 0.3*eV) < 1e-8

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

