def functional():
  from lada.vff import Functional
  vff = Functional()
  vff["In", "As"] = 2.62332, 21.6739, -112.0, 150.0
  vff["Ga", "As"] = 2.44795, 32.1530, -105.0, 150.0
  vff["As", "Ga", "As"] = "tet", -4.099, 9.3703
  vff["Ga", "As", "Ga"] = "tet", -4.099, 9.3703
  vff["In", "As", "In"] = "tet", -5.753, 5.7599
  vff["As", "In", "As"] = "tet", -5.753, 5.7599
  vff["Ga", "As", "In"] = -0.35016, -4.926, 7.5651

  return vff
def test_gradients(epsilon = 1e-4):
  from numpy import identity, abs, all, dot, array, sqrt
  from lada.crystal.binary import zinc_blende
  from quantities import eV, angstrom, newton, meter

  vff = functional()

  structure = zinc_blende()
  structure[0].type = 'In'
  structure[1].type = 'As'
  structure.scale = 6.5 #2.62332 * 2 / sqrt(3)  / 0.529177

  out = vff(structure)
  for atom, check in zip(structure, out):
    oldpos = atom.pos.copy()
    
    for dir in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
      dir = array(dir) / sqrt(dot(dir, dir))
      atom.pos = epsilon * dir + oldpos
      xplus = vff(structure).energy
      atom.pos = -epsilon * dir + oldpos
      xminus = vff(structure).energy

      assert abs(xplus - xminus) < 1e-8


  strain = array([[1e0, 0.1, 0], [0.1, 1e0, 0], [0, 0, 1e0]])
  structure.cell = dot(strain, structure.cell)
  for atom in structure: atom.pos = dot(strain, atom.pos)
  out = vff(structure)
  for atom, check in zip(structure, out):
    oldpos = atom.pos.copy()
    
    for dir in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
      dir = array(dir) / sqrt(dot(dir, dir))
      atom.pos = epsilon * dir + oldpos
      xplus = vff(structure).energy
      atom.pos = -epsilon * dir + oldpos
      xminus = vff(structure).energy

      deriv = (xplus - xminus) / (epsilon*angstrom)
      print deriv * 16.0217733 , dot(check.gradient, dir)


if __name__ == '__main__':
  test_gradients(1e-9)
