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

def test_inas():
  from numpy import identity, abs, all, dot, array
  from lada.crystal.binary import zinc_blende
  from quantities import eV

  vff = functional()

  structure = zinc_blende()
  structure[0].type = 'In'
  structure[1].type = 'As'
  structure.scale = 6.5 #2.62332 * 2 / sqrt(3)  / 0.529177

  out = vff(structure)
  assert abs(out.energy - 0.349587514524 * eV) < 1e-8
  assert all(abs(out.stress - identity(3) * 45.06322056) < 1e-8)
  assert all(abs(out[0].gradient) < 1e-8)
  assert all(abs(out[1].gradient) < 1e-8)

  epsilon = array([[1e0, 0.1, 0], [0.1, 1e0, 0], [0, 0, 1e0]])
  structure.cell = dot(epsilon, structure.cell)
  for atom in structure: atom.pos = dot(epsilon, atom.pos)
  out = vff(structure)
  assert abs(out.energy - 0.527010542895 * eV) < 1e-8
  assert all(abs(out.stress - [[  2.73218316e+01,   3.21556841e+01,  0],
                               [  3.21556841e+01,   2.73218316e+01,  0],
                               [ 0, 0,   2.01929521e+01]]) < 1e-6)
  assert all(abs(out[0].gradient - [0, 0, 2.27456489e+02]) < 1e-6)
  assert all(abs(out[1].gradient - [0, 0, -2.27456489e+02]) < 1e-6)

def test_ingaas():
  from numpy import abs, all, dot, array
  from quantities import eV
  from lada.crystal import Structure

  vff = functional()

  structure = Structure( 10.0, 0.5, 0.5, 
                         0.00, 0.0, 0.5,
                         0.00, 0.5, 0.0, scale=6.5 )\
                       .add_atom(pos=(0.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(0.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(1.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(1.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(2.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(2.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(3.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(3.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(4.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(4.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(5.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(5.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(6.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(6.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(7.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(7.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(8.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(8.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(9.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(9.25, 0.25, 0.25), type="As")  

  epsilon = array([[1e0, 0.1, 0], [0.1, 1e0, 0], [0, 0, 1e0]])
  structure.cell = dot(epsilon, structure.cell)
  for atom in structure: atom.pos = dot(epsilon, atom.pos)
  out = vff(structure)
  assert abs(out.energy - 12.7962077582*eV) < 1e-8
  assert all(abs(out.stress - [[ 767.82859534,  529.5647546 ,    2.75173374],
       [ 529.5647546 ,  767.82859534,   -2.75173374],
       [   2.75173374,   -2.75173374,  661.42954253]]) < 1e-6)
  gradients = array([u.gradient for u in out])
  check_gradients = array( 
      [[  1.56319402e-13,  -1.19726451e-12,   2.34055741e+02],
       [  0.00000000e+00,  -1.98241423e-12,  -2.34055741e+02],
       [  1.29119153e+01,  -1.29119153e+01,   2.09353996e+02],
       [ -1.21502808e-12,   4.83169060e-13,  -2.34055741e+02],
       [  1.59658154e+01,  -1.59658154e+01,   2.12736578e+02],
       [  7.96493276e+01,  -7.96493276e+01,  -3.44662597e+02],
       [ -3.95230085e+00,   3.95230085e+00,   2.41290573e+02],
       [ -1.77635684e-14,   1.90070182e-13,  -2.27456489e+02],
       [ -1.57183976e+01,   1.57183976e+01,   2.37666101e+02],
       [ -7.59444445e+01,   7.59444445e+01,  -1.19574166e+02],
       [  1.59658154e+01,  -1.59658154e+01,   2.12736578e+02],
       [  7.96493276e+01,  -7.96493276e+01,  -3.44662597e+02],
       [ -3.95230085e+00,   3.95230085e+00,   2.41290573e+02],
       [ -2.49222865e-12,  -4.83879603e-12,  -2.27456489e+02],
       [ -2.86303129e+01,   2.86303129e+01,   2.62367846e+02],
       [ -7.59444445e+01,   7.59444445e+01,  -1.19574166e+02],
       [ -1.70530257e-13,  -6.89226454e-13,   2.34055741e+02],
       [  1.72128978e-12,   2.13162821e-14,  -2.34055741e+02],
       [  7.10542736e-13,   1.53477231e-12,   2.34055741e+02],
       [ -1.42108547e-13,   1.54187774e-12,  -2.34055741e+02]] )
  assert all(abs(gradients - check_gradients) < 1e-6)

if __name__ == '__main__': 
  test_inas()
  test_ingaas()
