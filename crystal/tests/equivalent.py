""" Checks crystal method `lada.crystal.equivalent`. """
def scale(A, B):
  """ Check changes to scale. """
  from lada.crystal.cppwrappers import equivalent
  B = B.copy()
  assert equivalent(A, B)
  B.scale = 3.0;
  assert not equivalent(A, B)
  assert equivalent(A, B, scale=False)
  B.cell *= 0.5;
  for atom in B: atom.pos *= 0.5
  assert not equivalent(A, B)
  assert equivalent(A, B, scale=False)

def motif(A, B):
  """ Check changes in motif. """
  from numpy import dot
  from lada.crystal.cppwrappers import cellinvariants
  for op in cellinvariants(A):
    B = B.copy()
    for atom in B:
      atom.pos = dot(op[:3] * atom.pos)
      scale(A, B)

def basis(A, B):
  """ Adds rotation and translation of cartesian basis. """
  from numpy import dot, pi
  from lada.crystal.cppwrappers import cellinvariants
  from lada.math import Translation, Rotation, transform
  from random import random
  motif(A, transform(B, Rotation(0.5 * pi, [1,0,0])))
  motif(A, transform(B, Rotation(-pi, [1,0,0])))
  motif(A, transform(B, Rotation(-0.13*pi, [1,0,0])))
  motif(A, transform(B, Translation([0.25, 0.25, 0.25])))
  motif(A, transform(B, Rotation(random()*2*pi, [1, 0, 0]) \
                        + Translation([random()-0.5,random()-0.5,random()-0.5])))
def decoration(A, B):
  """ Adds changes to the motif. """
  
def test_fcc():
  """ Test fcc space-group and equivalents """
  from numpy import all, abs, identity, dot
  from numpy.linalg import inv, det
  from lada.crystal.cppwrappers import space_group, Structure, equivalent, transform
  from lada.math import is_integer

  structure = Structure([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]], m=True).add_atom(0,0,0,"Si", m=True)
  ops = space_group(structure)
  assert len(ops) == 48
  for op in ops:
    assert op.shape == (4, 3)
    assert all(abs(op[3, :]) < 1e-8) 

    other = transform(structure, op)
    assert all(abs(dot(op[:3], structure.cell)-other.cell) < 1e-8)
    assert getattr(other, 'm', False) 
    for a, atom in zip(structure, other):
      assert all(abs(dot(op[:3], a.pos) + op[3] - atom.pos) < 1e-8)
      assert a.type == atom.type
      assert getattr(atom, 'm', False) 

    assert equivalent(structure, other, cartesian=False)
    assert equivalent(other, structure, cartesian=False)

def test_b5(u):
  """ Test b5 space-group and equivalents """
  from numpy import all, abs, identity, dot
  from numpy.linalg import inv, det
  from lada.crystal.cppwrappers import space_group, Structure, equivalent, transform
  from lada.math import is_integer

  x, y = u, 0.25-u
  structure = Structure([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]) \
                       .add_atom(5.000000e-01, 5.000000e-01, 5.000000e-01, "A") \
                       .add_atom(5.000000e-01, 2.500000e-01, 2.500000e-01, "A") \
                       .add_atom(2.500000e-01, 5.000000e-01, 2.500000e-01, "A") \
                       .add_atom(2.500000e-01, 2.500000e-01, 5.000000e-01, "A") \
                       .add_atom(8.750000e-01, 8.750000e-01, 8.750000e-01, "B") \
                       .add_atom(1.250000e-01, 1.250000e-01, 1.250000e-01, "B") \
                       .add_atom(     x,     x,     x, "X") \
                       .add_atom(     x,     y,     y, "X") \
                       .add_atom(     y,     x,     y, "X") \
                       .add_atom(     y,     y,     x, "X") \
                       .add_atom(    -x,    -x,    -x, "X") \
                       .add_atom(    -x,    -y,    -y, "X") \
                       .add_atom(    -y,    -x,    -y, "X") \
                       .add_atom(    -y,    -y,    -x, "X") 
  ops = space_group(structure)
  assert len(ops) == 48
  for op in ops:
    assert op.shape == (4, 3)

    other = transform(structure, op)
    assert all(abs(dot(op[:3], structure.cell)-other.cell) < 1e-8)
    for a, atom in zip(structure, other):
      assert all(abs(dot(op[:3], a.pos) + op[3] - atom.pos) < 1e-8)
      assert a.type == atom.type

    assert equivalent(structure, other, cartesian=False)
    assert equivalent(other, structure, cartesian=False)
    
     

if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test_fcc()
  test_b5(0.25)
  test_b5(0.36)
