def functional():
  from lada.vff.vff import Vff
  vff = Vff()
  vff["In", "As"] = 2.62332, 21.6739, -112.0, 150.0
  vff["Ga", "As"] = 2.44795, 32.1530, -105.0, 150.0
  vff["As", "Ga", "As"] = "tet", -4.099, 9.3703
  vff["Ga", "As", "Ga"] = "tet", -4.099, 9.3703
  vff["In", "As", "In"] = "tet", -5.753, 5.7599
  vff["As", "In", "As"] = "tet", -5.753, 5.7599
  vff["Ga", "As", "In"] = -0.35016, -4.926, 7.5651

  return vff

def test_rep_vff():
  from numpy import all, abs
  from lada.vff.vff import Vff
  from lada.vff import exec_input

  a = Vff()
  input = exec_input(repr(a))
  assert len(input.vff._parameters) == 0

  a = functional()
  input = exec_input(repr(a))
  assert ("In", "As") in input.vff
  assert all(abs(input.vff["In", "As"] - a["In", "As"]) < 1e-8)
  assert all(abs(input.vff["In", "As", 'Ga'] - a["Ga", "As", 'In']) < 1e-8)

def test_pickle_vff():
  from pickle import loads, dumps
  from numpy import all, abs
  from lada.vff.vff import Vff

  a = Vff()
  b = loads(dumps(a)) 
  assert len(b._parameters) == 0

  a = functional()
  b = loads(dumps(a)) 
  assert ("In", "As") in b
  assert all(abs(b["In", "As"] - a["In", "As"]) < 1e-8)
  assert all(abs(b["In", "As", 'Ga'] - a["Ga", "As", 'In']) < 1e-8)

if __name__ == '__main__':
  test_rep_vff()
  test_pickle_vff()
