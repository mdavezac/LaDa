def test_random_cell():
  """ Tests whether scheduling jobs works on known failure cases. """
  from random import seed
  from numpy import abs, all, identity, array
  from lada.ga.xgsgo import random_cell

  assert all(abs(random_cell(cubic=True)-identity(3, dtype='float64'))<1e-8) 
  seed(1)
  assert all(abs(random_cell() - array([[-0.7074914 , -0.20268559, -0.13267866],
       [ 0.        ,  0.18178436,  1.27104087],
       [ 0.        ,  1.18001767,  0.        ]])) < 1e-8)

def test_random_structure():
  from random import seed
  from numpy import abs, all, identity, array
  from lada.ga.xgsgo import random_structure, populate_anion_cation, populate
  from lada.crystal import Structure

  result = Structure( -0.707491, -0.202686, -0.132679,\
                       0, 0.181784, 1.27104,\
                       0, 1.18002, 0,\
                       scale=4.34424 )\
              .add_atom(0, 0, 0, 'Si')\
              .add_atom(-0.576963, 0.648609, 0.674296, 'Ge')\
              .add_atom(-0.592216, 0.0779076, 0.505722, 'Si')\
              .add_atom(-0.0858174, 0.570701, 0.168574, 'Si')\
              .add_atom(-0.239546, 1.03773, 0.84287, 'Ge')\
              .add_atom(-0.528006, 0.933856, 0.168574, 'Ge')\
              .add_atom(-0.394825, 0.337392, 1.01144, 'Ge')\
              .add_atom(-0.380075, 0.233516, 0.337148, 'Si')
  seed(1)
  struc = random_structure(8)
  assert all(abs(struc.cell-result.cell) < 1e-4)
  for a, b in zip(result, struc):
    assert b.type == 'anion'
    assert all(abs(a.pos-b.pos) < 1e-4)

  populate_anion_cation(struc, {'Si': 4, 'Ge': 4}, ['Si'])
  for a, b in zip(result, struc):
    assert a.type == b.type
  
  populate(struc, {'Si': 4, 'Ge': 4})
  # shufle is not really that great it seems...
  for i, a in enumerate(struc):
    assert a.type == ('Si' if i < 4 else 'Ge')

if __name__ == "__main__":
  from sys import argv, path
  from os.path import abspath
  if len(argv) > 1: path.extend(argv[1:])

  test_random_cell()
  test_random_structure()
