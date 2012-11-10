def functional():
  from lada.vff.functional import Functional
  vff = Functional()
  vff["In", "As"] = 2.62332, 21.6739, -112.0, 150.0
  vff["Ga", "As"] = 2.44795, 32.1530, -105.0, 150.0
  vff["As", "Ga", "As"] = "tet", -4.099, 9.3703
  vff["Ga", "As", "Ga"] = "tet", -4.099, 9.3703
  vff["In", "As", "In"] = "tet", -5.753, 5.7599
  vff["As", "In", "As"] = "tet", -5.753, 5.7599
  vff["Ga", "As", "In"] = -0.35016, -4.926, 7.5651

  return vff

def test_inas(epsilon = 1e-4):
  from numpy import abs, dot, array, sqrt
  from lada.crystal.binary import zinc_blende
  from quantities import a0


  vff = functional()

  structure = zinc_blende()
  structure[0].type = 'In'
  structure[1].type = 'As'
  structure.scale = 2.62332 * 2 / sqrt(3)  * a0
  print structure.scale

  out = vff(structure)
  print out
      
if __name__ == '__main__':
  test_inas(1e-9)
